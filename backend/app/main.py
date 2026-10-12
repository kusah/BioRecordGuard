import os, hashlib
from io import StringIO
import requests
from Bio import SeqIO
from dotenv import load_dotenv
from fastapi import FastAPI, Depends, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy import create_engine, Column, Integer, String, Float, DateTime, JSON, desc
from sqlalchemy.orm import declarative_base, sessionmaker, Session
from sqlalchemy.sql import func

load_dotenv()
DB_URL = os.getenv('DATABASE_URL', 'sqlite:///./biorecordguard.db')
engine = create_engine(DB_URL, connect_args={'check_same_thread': False} if DB_URL.startswith('sqlite') else {}, pool_pre_ping=True)
SessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False)
Base = declarative_base()

class AuditRun(Base):
    __tablename__ = 'audit_runs'
    id = Column(Integer, primary_key=True, index=True)
    query_accession = Column(String(40), index=True, nullable=False)
    protein_name = Column(String(255))
    organism = Column(String(255))
    audit_score = Column(Float, nullable=False)
    status = Column(String(40), nullable=False)
    result_json = Column(JSON, nullable=False)
    created_at = Column(DateTime(timezone=True), server_default=func.now())

Base.metadata.create_all(bind=engine)
app = FastAPI(title='BioRecordGuard API', version='0.1.0')
app.add_middleware(CORSMiddleware, allow_origins=['http://localhost:5173','http://127.0.0.1:5173'], allow_credentials=True, allow_methods=['*'], allow_headers=['*'])

def get_db():
    db = SessionLocal()
    try: yield db
    finally: db.close()

def fetch_uniprot(accession):
    url = f'https://rest.uniprot.org/uniprotkb/{accession}.json'
    r = requests.get(url, timeout=20, headers={'User-Agent':'BioRecordGuard/0.1 student-demo'})
    r.raise_for_status(); d = r.json()
    desc = d.get('proteinDescription', {})
    name = desc.get('recommendedName', {}).get('fullName', {}).get('value')
    if not name:
        subs = desc.get('submissionNames', [])
        name = subs[0].get('fullName', {}).get('value') if subs else 'Not available'
    seq = d.get('sequence', {}).get('value', '')
    refs = [x.get('id') for x in d.get('uniProtKBCrossReferences', []) if x.get('database') == 'RefSeq' and x.get('id')]
    return {'source':'UniProtKB','accession':d.get('primaryAccession',accession),'name':name or 'Not available',
      'organism':d.get('organism',{}).get('scientificName','Not available'),'organism_id':d.get('organism',{}).get('taxonId'),
      'sequence':seq,'sequence_length':d.get('sequence',{}).get('length',len(seq)),
      'sequence_sha256':hashlib.sha256(seq.upper().encode()).hexdigest() if seq else None,
      'refseq_ids':refs,'source_url':url,'annotation_score':d.get('annotationScore'),'entry_type':d.get('entryType')}

def fetch_ncbi(refseq_id):
    params={'db':'protein','id':refseq_id,'rettype':'gp','retmode':'text','tool':'BioRecordGuard'}
    email=os.getenv('NCBI_EMAIL')
    if email: params['email']=email
    r=requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',params=params,timeout=20,headers={'User-Agent':'BioRecordGuard/0.1 student-demo'})
    r.raise_for_status(); raw=r.text.strip()
    if not raw or 'LOCUS' not in raw[:300]: raise ValueError('NCBI did not return a GenPept record')
    rec=SeqIO.read(StringIO(raw),'genbank'); seq=str(rec.seq)
    return {'source':'NCBI Protein','accession':rec.id,'name':rec.description or 'Not available',
      'organism':rec.annotations.get('organism','Not available'),'sequence':seq,'sequence_length':len(seq),
      'sequence_sha256':hashlib.sha256(seq.upper().encode()).hexdigest() if seq else None,
      'source_url':f'https://www.ncbi.nlm.nih.gov/protein/{rec.id}'}

def audit_records(uni,ncbi):
    checks=[]; findings=[]
    def add(label, applicable, passed, detail, weight=1):
        checks.append({'label':label,'applicable':applicable,'passed':passed,'detail':detail,'weight':weight})
        if applicable and not passed: findings.append({'severity':'review','rule':label,'detail':detail})
    add('Sequence length',ncbi is not None,ncbi is not None and uni['sequence_length']==ncbi['sequence_length'],
        f"UniProt length: {uni['sequence_length']}; NCBI length: {ncbi['sequence_length']}." if ncbi else 'No comparable NCBI sequence was retrieved.',2)
    exact=bool(ncbi and uni.get('sequence') and ncbi.get('sequence') and uni['sequence'].upper()==ncbi['sequence'].upper())
    add('Exact sequence match',ncbi is not None and bool(uni.get('sequence')) and bool(ncbi and ncbi.get('sequence')),exact,
        'Sequences are identical.' if exact else ('Sequences differ; this can reflect isoforms, versions, fragments or different record definitions.' if ncbi else 'No comparable NCBI sequence was retrieved.'),3)
    org_match=bool(ncbi and uni.get('organism','').casefold()==ncbi.get('organism','').casefold())
    add('Organism metadata',ncbi is not None,org_match,
        f"UniProt organism: {uni.get('organism')}; NCBI organism: {ncbi.get('organism')}." if ncbi else 'No comparable NCBI organism metadata.',2)
    name_ok=bool(uni.get('name') and uni['name']!='Not available')
    add('Protein name present',True,name_ok,'Protein name is present.' if name_ok else 'Protein name is missing.',1)
    seq_ok=bool(uni.get('sequence'))
    add('UniProt sequence present',True,seq_ok,'Sequence is present.' if seq_ok else 'Sequence is missing.',2)
    active=[c for c in checks if c['applicable']]; total=sum(c['weight'] for c in active)
    score=round(100*sum(c['weight'] for c in active if c['passed'])/total) if total else 0
    return {'score':score,'status':'No review flags' if not findings else 'Review recommended',
      'score_note':"Internal rule-based audit score; not an official database quality score or probability of correctness.",
      'checks':checks,'findings':findings}

@app.get('/api/health')
def health(): return {'status':'ok','app':'BioRecordGuard'}

@app.post('/api/audit/{accession}')
def run_audit(accession: str, db: Session=Depends(get_db)):
    accession=accession.strip()
    if not accession or len(accession)>40 or any(c.isspace() for c in accession): raise HTTPException(400,'Enter a valid accession.')
    try: uni=fetch_uniprot(accession)
    except requests.HTTPError as e:
        code=e.response.status_code if e.response is not None else 502
        raise HTTPException(404 if code==404 else 502, f'UniProt request failed (HTTP {code}). Check accession and connection.')
    except requests.RequestException as e: raise HTTPException(502,f'Could not reach UniProt: {e}')
    ncbi=None; message='No RefSeq cross-reference was available in this UniProt record.'
    for ref in uni.get('refseq_ids',[])[:3]:
        try: ncbi=fetch_ncbi(ref); message=f'Retrieved linked NCBI Protein record {ncbi["accession"]}.'; break
        except (requests.RequestException,ValueError) as e: message=f'Linked RefSeq {ref} could not be retrieved: {e}'
    audit=audit_records(uni,ncbi)
    result={'query_accession':accession,'uniprot':{k:v for k,v in uni.items() if k not in ('sequence','refseq_ids')},
      'ncbi':({k:v for k,v in ncbi.items() if k!='sequence'} if ncbi else None),'ncbi_message':message,'audit':audit}
    row=AuditRun(query_accession=accession,protein_name=uni.get('name'),organism=uni.get('organism'),audit_score=audit['score'],status=audit['status'],result_json=result)
    db.add(row); db.commit(); db.refresh(row)
    result['audit_id']=row.id; result['created_at']=row.created_at.isoformat() if row.created_at else None
    return result

@app.get('/api/history')
def history(limit:int=10,db:Session=Depends(get_db)):
    rows=db.query(AuditRun).order_by(desc(AuditRun.id)).limit(max(1,min(limit,50))).all()
    return [{'id':r.id,'query_accession':r.query_accession,'protein_name':r.protein_name,'organism':r.organism,'audit_score':r.audit_score,'status':r.status,'created_at':r.created_at.isoformat() if r.created_at else None} for r in rows]

@app.get('/api/history/{audit_id}')
def history_item(audit_id:int,db:Session=Depends(get_db)):
    row=db.query(AuditRun).filter(AuditRun.id==audit_id).first()
    if not row: raise HTTPException(404,'Audit run not found.')
    result=dict(row.result_json); result['audit_id']=row.id; result['created_at']=row.created_at.isoformat() if row.created_at else None
    return result
