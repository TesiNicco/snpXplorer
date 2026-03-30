import os, sqlite3, re
import pandas as pd
from typing import List, Tuple, Dict

DIST_MAX = 1_000_000   # only consider LD within this distance
R2_SCALE = 1000        # r² is stored as int(r²*1000)
MMAP_BYTES = 1 << 30

def _open_db(db_path: str) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")
    conn.execute(f"PRAGMA mmap_size={MMAP_BYTES};")
    return conn

def _detect_chr_column(conn: sqlite3.Connection) -> str:
    cols = {r[1].lower() for r in conn.execute("PRAGMA table_info(variants);")}
    for cand in ("chr","chrom","chromosome"):
        if cand in cols: return cand
    return None

def _infer_chr_from_path(path: str) -> str:
    m = re.search(r"(chr[\w\-]+)", os.path.basename(path).lower())
    return m.group(1) if m else "chr"

def get_variants(conn: sqlite3.Connection, chr_name: str) -> List[Tuple[int,int,str]]:
    """Return list of (id,pos,uniq) for variants in this chr ordered by pos."""
    chr_col = _detect_chr_column(conn)
    if chr_col:
        q = f"SELECT id,pos,uniq FROM variants WHERE {chr_col}=? ORDER BY pos;"
        rows = conn.execute(q,(chr_name,)).fetchall()
    else:
        q = "SELECT id,pos,uniq FROM variants ORDER BY pos;"
        rows = conn.execute(q).fetchall()
    return rows

def get_pairwise_ld(conn: sqlite3.Connection, v_ids: List[int], r2m_min: int) -> Dict[Tuple[int,int], float]:
    """Fetch LD for pairs in v_ids (within DIST_MAX). Return {(id1,id2): r2} dict."""
    idset = set(v_ids)
    out = {}
    q = f"""
    SELECT v1,v2,r2_milli
    FROM ld
    WHERE r2_milli >= ? AND dist_bp <= ?
    """
    for v1,v2,r2 in conn.execute(q,(r2m_min,DIST_MAX)):
        if v1 in idset and v2 in idset:
            out[(min(v1,v2),max(v1,v2))] = r2/ R2_SCALE
    return out

def write_haplotype_blocks(ld_db: str, out_folder: str="../Data/databases/HaploBlocks", r2_min: float=0.8) -> str:
    os.makedirs(out_folder, exist_ok=True)
    chr_name = _infer_chr_from_path(ld_db)
    out_name = f"{chr_name}_haploblocks_r2_{str(r2_min).replace('.','p')}.sqlite"
    out_path = os.path.join(out_folder,out_name)
    conn_in = _open_db(ld_db)
    conn_out = sqlite3.connect(out_path)
    conn_out.execute("PRAGMA journal_mode=WAL;")
    conn_out.execute("PRAGMA synchronous=NORMAL;")
    conn_out.execute("""
    CREATE TABLE IF NOT EXISTS haplotypes(
      chr TEXT,
      haplo_id INTEGER,
      start_pos INTEGER,
      end_pos INTEGER,
      snps TEXT,
      PRIMARY KEY(chr,haplo_id)
    );""")
    r2m_min = int(r2_min*R2_SCALE)
    variants = get_variants(conn_in, chr_name)
    print(f"[{chr_name}] {len(variants)} variants")
    haplo_id, blocks = 1, []
    i = 0
    while i < len(variants):
        block = [variants[i]]
        j = i+1
        while j < len(variants):
            cand_block = block + [variants[j]]
            ids = [v[0] for v in cand_block]
            ldmap = get_pairwise_ld(conn_in, ids, r2m_min)
            # check all pairs
            ok = True
            for a in range(len(cand_block)):
                for b in range(a+1,len(cand_block)):
                    id1,id2 = cand_block[a][0], cand_block[b][0]
                    if (min(id1,id2),max(id1,id2)) not in ldmap:
                        ok = False
                        break
                if not ok: break
            if not ok: break
            block = cand_block
            j += 1
        poses = [v[1] for v in block]
        uniqs = [v[2] for v in block]
        blocks.append((chr_name,haplo_id,min(poses),max(poses),",".join(uniqs)))
        haplo_id += 1
        i = j  # start next block after last included SNP
    conn_out.executemany("INSERT INTO haplotypes(chr,haplo_id,start_pos,end_pos,snps) VALUES (?,?,?,?,?)",blocks)
    conn_out.commit()
    conn_in.close(); conn_out.close()
    print(f"[{chr_name}] wrote {len(blocks)} blocks to {out_path}")
    return out_path

def haplotypes_in_interval(haplo_db: str, chr_name: str, start: int, end: int):
    conn = sqlite3.connect(haplo_db)
    q = """SELECT haplo_id,start_pos,end_pos,snps FROM haplotypes
           WHERE chr=? AND end_pos>=? AND start_pos<=? ORDER BY haplo_id;"""
    rows = conn.execute(q,(chr_name,start,end)).fetchall()
    conn.close()
    return pd.DataFrame(rows,columns=["haplo_id","start_pos","end_pos","snps"])

out_db = write_haplotype_blocks("../Data/databases/LD_db/ld_chr16.sqlite", r2_min=0.8)
df = haplotypes_in_interval(out_db, "chr16", 81563046, 82246253)
print(df)