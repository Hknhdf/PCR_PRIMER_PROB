"""
conservation.py
Coklu dizi hizalamasi ile korunmus bolge (conservation) analizi.

Yontem:
  1. NCBI'dan ayni organizmaya ait aday diziler cekilir.
  2. Adaylar 15-mer benzerligi ile on filtrelenir (homolog olmayan kayitlar elenir).
  3. Referans bolge (hedef bolge + kanatlar) alinir; her dizide 24-mer capa aranir.
     Capa bulunursa yalnizca o pencere global hizalanir (hizli), bulunamazsa lokal hizalama.
  4. Hizalama bloklari referans koordinatlarina izdusurulur.
  5. Her pozisyon icin "ayni bazda olan dizi orani" hesaplanir; veri yoksa None.

Harici arac (MAFFT/MUSCLE) gerektirmez.
"""

from __future__ import annotations

from Bio import Entrez, SeqIO
from Bio.Align import PairwiseAligner

ANCHOR_K = 24
FLANK = 120
MIN_IDENTITY = 0.70
MIN_COVERAGE = 0.50
KMER_K = 15
KMER_STEP = 3
MIN_KMER_SIMILARITY = 0.08


def _make_aligner(mode="global"):
    al = PairwiseAligner()
    al.mode = mode
    al.match_score = 2.0
    al.mismatch_score = -1.0
    al.open_gap_score = -6.0
    al.extend_gap_score = -1.0
    return al


def _kmers(seq, k=KMER_K, step=KMER_STEP):
    return {seq[i:i + k] for i in range(0, max(0, len(seq) - k + 1), step)}


def kmer_similarity(reference: str, subject: str, k=KMER_K, step=KMER_STEP) -> float:
    """Referansin k-mer'larinin ne kadari subject'te var (0..1).
    Buyuk genomlarda set kurmadan tarama yapilir (bellek dostu)."""
    ref = str(reference or "").upper()
    sub = str(subject or "").upper()
    if len(ref) < k or len(sub) < k:
        return 0.0
    ref_kmers = _kmers(ref, k, step)
    if not ref_kmers:
        return 0.0
    if len(sub) > 3 * len(ref):
        # Büyük genom: kmer seti kurmak yerine referans kmer'larını genomda ara.
        # (Adım/frame kayması sorununu önler.)
        if len(ref_kmers) <= 20_000:
            hits = sum(1 for km in ref_kmers if sub.find(km) >= 0)
            return hits / len(ref_kmers)
        found = set()
        for frame in range(3):
            for i in range(frame, len(sub) - k + 1, step):
                km = sub[i:i + k]
                if km in ref_kmers:
                    found.add(km)
        return len(found) / len(ref_kmers)
    return len(ref_kmers & _kmers(sub, k, step)) / len(ref_kmers)


def region_keywords(record, start, end, limit=3):
    """
    Hedef bölgeyle örtüşen gen/CDS ürün adlarını döndürür (NCBI arama anahtarı olarak).
    ör. 'nucleoprotein', 'envelope', 'NS5'
    """
    keywords = []
    try:
        for feature in record.features:
            if feature.type not in ("CDS", "gene", "mRNA"):
                continue
            f_start = int(feature.location.start)
            f_end = int(feature.location.end)
            if not (f_start <= end and f_end >= start):
                continue
            for key in ("gene", "product"):
                for value in feature.qualifiers.get(key, []):
                    value = str(value).strip()
                    if value and len(value) <= 40 and value not in keywords:
                        keywords.append(value)
    except Exception:
        return []
    return keywords[:limit]


def fetch_related_sequences(
    organism,
    exclude_accessions=(),
    max_n=12,
    reference_seq=None,
    min_similarity=MIN_KMER_SIMILARITY,
    keywords=(),
    max_candidates=None,
    min_subject_len=None,
):
    """
    Ayni organizmaya ait nukleotid kayitlarini ceker, homolog olanlari dondurur.
    keywords: hedef gendeki urun adlari (arama hassasiyetini artirir).
    reference_seq: k-mer benzerligi filtresi icin referans dizi.
    max_candidates: buyuk genomlarda indirilecek kayit sayisini sinirlar.
    Donus: [(accession, description, sequence_str), ...]
    """
    if not organism:
        return []

    excluded = {str(a).upper() for a in exclude_accessions}
    ref = str(reference_seq).upper() if reference_seq else None
    cap = int(max_candidates) if max_candidates else max(60, max_n * 6)
    # Buyuk genomlarda istek basina kayit sayisi kucultulur (yanit boyutu / bellek)
    chunk_size = 3 if max_candidates else 50

    def _search(term):
        try:
            with Entrez.esearch(db="nucleotide", term=term, retmax=cap) as handle:
                found = Entrez.read(handle)
            return [int(i) for i in found.get("IdList", [])]
        except Exception:
            return []

    def _collect(ids):
        picked = []
        for start in range(0, len(ids), chunk_size):
            chunk = ids[start:start + chunk_size]
            try:
                with Entrez.efetch(
                    db="nucleotide", id=",".join(str(i) for i in chunk),
                    rettype="fasta", retmode="text",
                ) as handle:
                    records = list(SeqIO.parse(handle, "fasta"))
            except Exception:
                continue

            for rec in records:
                acc = rec.id.split(".")[0].upper()
                if acc in excluded:
                    continue
                seq = str(rec.seq).upper()
                if len(seq) < 60 or set(seq) - set("ACGTN"):
                    continue
                if min_subject_len and len(seq) < min_subject_len:
                    continue
                if ref:
                    sim = kmer_similarity(ref, seq)
                    if sim < min_similarity:
                        continue
                else:
                    sim = 1.0
                picked.append((sim, rec.id, rec.description, seq))
                if len(picked) >= max_n:
                    break
            if len(picked) >= max_n:
                break

        picked.sort(key=lambda x: x[0], reverse=True)
        return [(cid, desc, seq) for _, cid, desc, seq in picked[:max_n]]

    # 1) Hedef gen anahtar kelimeleriyle ara (daha isabetli)
    if keywords:
        kw_terms = " OR ".join(f'"{k}"' for k in keywords[:3])
        hits = _collect(_search(f'"{organism}"[Organism] AND ({kw_terms})'))
        if hits:
            return hits

    # 2) Organizma genelinde ara (yedek)
    return _collect(_search(f'"{organism}"[Organism]'))


def _find_anchor(ref_slice: str, subject: str, k: int = ANCHOR_K):
    """ref_slice'i subject icinde kabaca konumlandir (offset: ref 0 -> subject)."""
    L = len(ref_slice)
    for frac in (0.5, 0.35, 0.65, 0.2, 0.8, 0.45, 0.6, 0.3):
        p = max(0, min(L - k, int(L * frac) - k // 2))
        probe = ref_slice[p:p + k]
        if len(probe) < k:
            continue
        idx = subject.find(probe)
        if idx >= 0:
            return idx - p
    return None


def _project_blocks(alignment):
    """Hizalama bloklarini (ref_index, subject_index) ciftlerine cevirir."""
    ref_blocks, subj_blocks = alignment.aligned
    pairs = []
    for (r0, r1), (s0, s1) in zip(ref_blocks, subj_blocks):
        length = min(int(r1 - r0), int(s1 - s0))
        for k in range(length):
            pairs.append((int(r0) + k, int(s0) + k))
    return pairs


def compute_conservation(
    ref_seq: str,
    region_start: int,
    region_end: int,
    subjects,
    flank: int = FLANK,
    min_identity: float = MIN_IDENTITY,
    min_coverage: float = MIN_COVERAGE,
):
    """
    Referans [region_start, region_end) bolgesi icin pozisyon bazli korunmusluk.

    Donus:
      {
        "conservation": [0..1 veya None] * (region_end-region_start),
        "coverage": [0..1] ayni uzunlukta,
        "mean": float (bilinen pozisyonlar uzerinden),
        "known": int (veri olan pozisyon sayisi),
        "n_used": int,
        "details": [{"id","description","identity","coverage"}...],
      }
    """
    ref_seq = str(ref_seq).upper()
    region_start = max(0, int(region_start))
    region_end = min(len(ref_seq), int(region_end))
    n = region_end - region_start
    if n <= 0:
        return {"conservation": [], "coverage": [], "mean": 0.0, "known": 0,
                "n_used": 0, "details": []}

    slice_start = max(0, region_start - flank)
    slice_end = min(len(ref_seq), region_end + flank)
    ref_slice = ref_seq[slice_start:slice_end]
    L = len(ref_slice)

    cover = [0] * L
    match = [0] * L
    details = []

    aligner_global = _make_aligner("global")
    aligner_local = _make_aligner("local")

    for item in subjects:
        try:
            seq_id, description, subject = item[0], item[1], item[2]
        except Exception:
            continue
        subject = str(subject).upper()
        if len(subject) < 60:
            continue

        offset = _find_anchor(ref_slice, subject)
        try:
            if offset is not None:
                win_start = max(0, offset - 250)
                win_end = min(len(subject), offset + L + 250)
                window = subject[win_start:win_end]
                alignments = aligner_global.align(ref_slice, window)
                if len(alignments) == 0:
                    continue
                pairs = _project_blocks(alignments[0])
                base_seq, base_shift = window, win_start
            else:
                if len(subject) > 40000:
                    continue
                alignments = aligner_local.align(ref_slice, subject)
                if len(alignments) == 0:
                    continue
                pairs = _project_blocks(alignments[0])
                base_seq, base_shift = subject, 0
        except Exception:
            continue

        if not pairs:
            continue

        # Bu dizi icin gecici sayim (kabul edilmezse ana sayaca eklenmez)
        sub_cover = [0] * L
        sub_match = [0] * L
        for r_idx, s_idx in pairs:
            if not (0 <= r_idx < L):
                continue
            s_abs = s_idx
            if not (0 <= s_abs < len(base_seq)):
                continue
            sub_cover[r_idx] += 1
            if base_seq[s_abs] == ref_slice[r_idx]:
                sub_match[r_idx] += 1

        counted = sum(1 for v in sub_cover if v)
        if counted == 0:
            continue
        identity = sum(sub_match) / counted
        coverage = counted / L
        if identity < min_identity or coverage < min_coverage:
            continue

        for i in range(L):
            if sub_cover[i]:
                cover[i] += 1
                match[i] += sub_match[i]

        details.append({
            "id": seq_id,
            "description": description,
            "identity": round(identity * 100, 1),
            "coverage": round(coverage * 100, 1),
        })

    idx0 = region_start - slice_start
    idx1 = idx0 + n
    conservation, cov_out = [], []
    known_vals = []
    for i in range(idx0, idx1):
        c = cover[i]
        if c == 0:
            conservation.append(None)
            cov_out.append(0.0)
        else:
            val = round(match[i] / c, 3)
            conservation.append(val)
            cov_out.append(round(c / len(details), 3) if details else 0.0)
            known_vals.append(val)

    mean = round(sum(known_vals) / len(known_vals), 3) if known_vals else 0.0
    return {
        "conservation": conservation,
        "coverage": cov_out,
        "mean": mean,
        "known": len(known_vals),
        "n_used": len(details),
        "details": details,
    }


def window_means(conservation, window=30):
    """Grafik icin pencere ortalamalari (veri olmayan pencereler None)."""
    vals = list(conservation or [])
    if not vals:
        return []
    out = []
    step = max(1, window // 3)
    for i in range(0, len(vals), step):
        chunk = [v for v in vals[i:i + window] if v is not None]
        out.append(round(sum(chunk) / len(chunk), 3) if chunk else None)
    return out


def best_conserved_windows(conservation, window=60, top_n=5, min_mean=0.85):
    """Ortalamasi en yuksek korunmus pencereler: [(start, end, mean)]"""
    vals = list(conservation or [])
    if len(vals) < window:
        return []

    scored = []
    for i in range(0, len(vals) - window + 1):
        chunk = [v for v in vals[i:i + window] if v is not None]
        if len(chunk) < window * 0.5:
            continue
        scored.append((i, i + window, sum(chunk) / len(chunk)))

    scored.sort(key=lambda x: x[2], reverse=True)
    picked = []
    for s, e, m in scored:
        if m < min_mean:
            break
        if any(not (e <= ps or s >= pe) for ps, pe, _ in picked):
            continue
        picked.append((s, e, round(m, 3)))
        if len(picked) >= top_n:
            break
    return picked
