# BioRecordGuard demo starter

Student prototype: retrieves a UniProtKB record, attempts to retrieve linked NCBI Protein via a RefSeq cross-reference, compares selected fields, stores audit runs, and displays a React dashboard.

## Backend (PowerShell)
```powershell
cd backend
python -m venv venv
.\venv\Scripts\Activate.ps1
pip install -r requirements.txt
copy .env.example .env
uvicorn app.main:app --reload
```
API docs: http://127.0.0.1:8000/docs

## Frontend (second terminal)
```powershell
cd frontend
npm install
npm run dev
```
Open http://localhost:5173. Demo accession: `P04637`.

## MySQL option
The default is SQLite for a quick, reliable demo. To use MySQL, create the database:
```sql
CREATE DATABASE biorecordguard CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
```
Then edit `backend/.env` with:
`DATABASE_URL=mysql+pymysql://YOUR_USER:YOUR_PASSWORD@localhost:3306/biorecordguard`
Add your email to `NCBI_EMAIL`, restart backend. Do not commit `.env`.

The audit score is an internal rule-based score, not an official quality score or probability. Differences may be legitimate due to isoforms, versions, fragments or record definitions. If a linked NCBI record is unavailable, the UI says so rather than inventing a match.
