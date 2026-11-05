"""
specificity.py
Primer/prob ozgulluk motoru:
  - seed-and-extend ile <=N mismatch toleransli baglanma bolgeleri
  - in-silico PCR (amplikon + prob konumu)
  - panel taramasi ve risk skoru
  - istege bagli NCBI BLAST URL API istemcisi
"""

from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from typing import Iterable

import requests
from Bio.Blast import NCBIXML

from pcr_engine import revcomp


@dataclass
class BindingSite:
    position: int                 # 0-based, sense strand
    mismatches: list[int]         # oligo pozisyonlari (0-based)
    n_mismatches: int
    three_prime_mm: int
    five_prime_mm: int
    strand: str = "+"

    def as_dict(self):
        return {
            "position": self.position,
            "position_1based": self.position + 1,
            "mismatches": self.mismatches,
            "n_mismatches": self.n_mismatches,
            "three_prime_mm": self.three_prime_mm,
            "five_prime_mm": self.five_prime_mm,
            "strand": self.strand,
        }


@dataclass
class PcrHit:
    forward: BindingSite
    reverse: BindingSite
    amplicon_start: int           # 0-based, inclusive
    amplicon_end: int             # 0-based, exclusive
    amplicon: str
    probe: BindingSite | None = None
    probe_offset: int | None = None
    extra: dict = field(default_factory=dict)

    @property
    def size(self):
        return self.amplicon_end - self.amplicon_start

    @property
    def total_mismatches(self):
        n = self.forward.n_mismatches + self.reverse.n_mismatches
        if self.probe is not None:
            n += self.probe.n_mismatches
        return n

    @property
    def three_prime_mismatches(self):
        return self.forward.three_prime_mm + self.reverse.three_prime_mm

    def as_dict(self):
        return {
            "forward": self.forward.as_dict(),
            "reverse": self.reverse.as_dict(),
            "amplicon_start": self.amplicon_start,
            "amplicon_end": self.amplicon_end,
            "amplicon_start_1based": self.amplicon_start + 1,
            "amplicon_end_1based": self.amplicon_end,
            "size": self.size,
            "probe": self.probe.as_dict() if self.probe else None,
            "probe_offset": self.probe_offset,
            "total_mismatches": self.total_mismatches,
            "three_prime_mismatches": self.three_prime_mismatches,
        }


def _seed_blocks(oligo: str, max_mismatches: int) -> list[tuple[int, str]]:
    """En az bir blogun mismatchsiz olmasini garanti eden bloklar."""
    L = len(oligo)
    n_blocks = min(L, max(1, int(max_mismatches) + 1))
    size = max(1, L // n_blocks)
    blocks = []
    start = 0
    while start < L:
        block = oligo[start:start + size]
        if block:
            blocks.append((start, block))
        start += size
    # cok kisa son bloklari birlestir
    if len(blocks) > 1 and len(blocks[-1][1]) < 4:
        off, blk = blocks.pop()
        prev_off, prev_blk = blocks[-1]
        blocks[-1] = (prev_off, prev_blk + blk)
    return blocks


def find_binding_sites(oligo: str, seq: str, max_mismatches: int = 3,
                       three_prime_side: str = "right") -> list[BindingSite]:
    """
    Oligonun sense dizideki baglanma bolgelerini bulur.
    three_prime_side:
      'right' -> oligonun 3' ucu pencerenin saginda (forward primer / prob)
      'left'  -> oligonun 3' ucu pencerenin solunda (reverse primerin revcomp hedefi)
    """
    oligo = str(oligo).upper()
    seq = str(seq).upper()
    L = len(oligo)
    if L == 0 or len(seq) < L:
        return []

    positions: dict[int, BindingSite] = {}
    for offset, block in _seed_blocks(oligo, max_mismatches):
        start = 0
        while True:
            idx = seq.find(block, start)
            if idx == -1:
                break
            pos = idx - offset
            if 0 <= pos <= len(seq) - L:
                window = seq[pos:pos + L]
                mm = [i for i, (a, b) in enumerate(zip(window, oligo)) if a != b]
                if len(mm) <= max_mismatches and pos not in positions:
                    if three_prime_side == "left":
                        three_mm = sum(1 for i in mm if i < 5)
                        five_mm = sum(1 for i in mm if i >= L - 5)
                    else:
                        three_mm = sum(1 for i in mm if i >= L - 5)
                        five_mm = sum(1 for i in mm if i < 5)
                    positions[pos] = BindingSite(
                        position=pos, mismatches=mm, n_mismatches=len(mm),
                        three_prime_mm=three_mm, five_prime_mm=five_mm,
                    )
            start = idx + 1
    return sorted(positions.values(), key=lambda s: (s.position, s.n_mismatches))


def in_silico_pcr(forward: str, reverse: str, seq: str, probe: str | None = None,
                  max_mismatches: int = 3, min_amp: int = 40,
                  max_amp: int = 1000) -> list[PcrHit]:
    """Tek bir dizide toleransli in-silico PCR."""
    forward = str(forward).upper()
    reverse = str(reverse).upper()
    seq = str(seq).upper()
    rev_target = revcomp(reverse)

    f_sites = find_binding_sites(forward, seq, max_mismatches, "right")
    r_sites = find_binding_sites(rev_target, seq, max_mismatches, "left")
    if not f_sites or not r_sites:
        return []

    hits: list[PcrHit] = []
    for f in f_sites:
        for r in r_sites:
            if r.position < f.position + len(forward):
                continue
            amp_start = f.position
            amp_end = r.position + len(reverse)
            size = amp_end - amp_start
            if size < min_amp or size > max_amp:
                continue
            amplicon = seq[amp_start:amp_end]
            probe_site = None
            probe_offset = None
            if probe:
                probe = str(probe).upper()
                for site in find_binding_sites(probe, amplicon, max_mismatches, "right"):
                    p_start = site.position
                    p_end = p_start + len(probe)
                    if p_start >= len(forward) and p_end <= len(amplicon) - len(reverse):
                        probe_site = site
                        probe_offset = p_start
                        break
                if probe_site is None:
                    # prob baglanmiyorsa yine de primer ciftini raporla
                    pass
            hits.append(PcrHit(
                forward=f, reverse=r,
                amplicon_start=amp_start, amplicon_end=amp_end,
                amplicon=amplicon, probe=probe_site, probe_offset=probe_offset,
            ))
    hits.sort(key=lambda h: (h.total_mismatches, h.three_prime_mismatches, h.size))
    return hits


def scan_sequences(forward: str, reverse: str, probe: str | None,
                   sequences: dict[str, str], max_mismatches: int = 3,
                   min_amp: int = 40, max_amp: int = 1000) -> list[dict]:
    """
    sequences: {isim/accession: dizi}
    Donus: her dizi icin en iyi hit + tum hit'ler.
    """
    records = []
    for name, seq in sequences.items():
        hits = in_silico_pcr(forward, reverse, seq, probe,
                             max_mismatches=max_mismatches,
                             min_amp=min_amp, max_amp=max_amp)
        records.append({
            "name": name,
            "length": len(seq),
            "n_hits": len(hits),
            "best": hits[0] if hits else None,
            "hits": hits,
        })
    return records


def risk_level(hit: PcrHit | None, is_target: bool = False) -> str:
    if hit is None:
        return "yok"
    if is_target:
        return "hedef"
    if hit.total_mismatches == 0:
        return "yuksek"
    if hit.total_mismatches <= 1 or (hit.total_mismatches <= 2 and hit.three_prime_mismatches == 0):
        return "orta"
    return "dusuk"


def specificity_score(target_hit: PcrHit | None, off_targets: Iterable[dict]) -> tuple[float, list[str]]:
    """
    0-100 ozgulluk skoru. Hedef disi amplikonlar ceza alir.
    off_targets: scan_sequences cikisi (is_target=False olanlar).
    """
    score = 100.0
    notes = []
    if target_hit is None:
        score = 0.0
        notes.append("Hedef dizide amplikon bulunamadi.")
    else:
        if target_hit.three_prime_mismatches:
            score -= 10
            notes.append(f"Hedefte 3' uc mismatch: {target_hit.three_prime_mismatches}")
        if target_hit.probe is None:
            score -= 10
            notes.append("Prob hedef amplikona baglanmiyor.")
    for rec in off_targets:
        best = rec.get("best")
        if best is None:
            continue
        risk = risk_level(best, is_target=False)
        if risk == "yuksek":
            score -= 40
            notes.append(f"{rec['name']}: mukemmel off-target amplikon ({best.size} bp)")
        elif risk == "orta":
            score -= 18
            notes.append(f"{rec['name']}: orta risk ({best.total_mismatches} mismatch, {best.size} bp)")
        elif risk == "dusuk":
            score -= 6
            notes.append(f"{rec['name']}: dusuk risk ({best.total_mismatches} mismatch)")
    return round(max(0.0, min(100.0, score)), 1), notes


# ---------------------------------------------------------------------------
# NCBI BLAST URL API (istege bagli; ag/gecikme durumunda hata dondurur)
# ---------------------------------------------------------------------------

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"


def _parse_blast_xml(xml_text: str, max_hits: int = 25) -> list[dict]:
    import io

    hits = []
    try:
        records = NCBIXML.parse(io.StringIO(xml_text))
        for record in records:
            for alignment in record.alignments[:max_hits]:
                for hsp in alignment.hsps[:3]:
                    hits.append({
                        "accession": alignment.accession,
                        "title": alignment.title[:160],
                        "length": alignment.length,
                        "evalue": hsp.expect,
                        "score": hsp.score,
                        "identity": round(100.0 * hsp.identities / max(1, hsp.align_length), 1),
                        "align_length": hsp.align_length,
                        "mismatches": int(hsp.align_length - hsp.identities),
                        "gaps": hsp.gaps,
                        "query_start": hsp.query_start,
                        "query_end": hsp.query_end,
                        "subject_start": hsp.sbjct_start,
                        "subject_end": hsp.sbjct_end,
                    })
    except Exception:
        return []
    hits.sort(key=lambda h: (h["evalue"], -h["identity"]))
    return hits[:max_hits]


def blast_oligo(oligo: str, database: str = "nt", max_wait: int = 150,
                human_only: bool = False, max_hits: int = 25) -> dict:
    """
    Tek oligoyu NCBI BLAST ile tarar. Basarisizlikta {'error': ...} dondurur.
    Not: resmi Primer-BLAST API'si yoktur; bu, blastn-short tabanli bir kontroludur.
    """
    oligo = str(oligo).upper()
    if len(oligo) < 15:
        return {"error": "BLAST icin oligo en az 15 nt olmali."}
    params = {
        "CMD": "Put",
        "PROGRAM": "blastn",
        "DATABASE": database,
        "QUERY": oligo,
        "TASK": "blastn-short",
        "WORD_SIZE": 7,
        "DUST": "no",
        "EXPECT": 1000,
        "HITLIST_SIZE": max(50, max_hits),
        "FORMAT_TYPE": "XML",
    }
    if human_only:
        params["ENTREZ_QUERY"] = "txid9606[ORGN]"
    try:
        put = requests.post(BLAST_URL, data=params, timeout=30)
        put.raise_for_status()
    except Exception as exc:
        return {"error": f"BLAST istegi gonderilemedi: {exc}"}

    m = re.search(r"RID = (\S+)", put.text)
    if not m:
        return {"error": "BLAST RID alinamadi (NCBI yogunlugu veya engel)."}
    rid = m.group(1)

    elapsed = 0
    while elapsed < max_wait:
        time.sleep(5)
        elapsed += 5
        try:
            poll = requests.get(BLAST_URL, params={
                "CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo",
            }, timeout=30)
            poll.raise_for_status()
        except Exception:
            continue
        if "Status=WAITING" in poll.text:
            continue
        if "Status=FAILED" in poll.text:
            return {"error": "BLAST islemi basarisiz oldu."}
        if "Status=UNKNOWN" in poll.text:
            return {"error": "BLAST RID bulunamadi."}
        if "Status=READY" in poll.text:
            try:
                res = requests.get(BLAST_URL, params={
                    "CMD": "Get", "RID": rid, "FORMAT_TYPE": "XML",
                    "FORMAT_OBJECT": "Alignment",
                }, timeout=60)
                res.raise_for_status()
            except Exception as exc:
                return {"error": f"BLAST sonucu alinamadi: {exc}"}
            return {"rid": rid, "hits": _parse_blast_xml(res.text, max_hits)}
    return {"error": f"BLAST {max_wait} s icinde tamamlanmadi (RID {rid})."}
