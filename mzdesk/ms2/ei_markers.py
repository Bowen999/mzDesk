"""Marker ions of GC-MS silyl derivatives in a 70 eV EI spectrum (step 3 of the EI MS Annotation workflow). Pure functions, no engine.

What the rule reads (all intensities are % of the base peak of the spectrum, ions found within +-TOLERANCE of the integer m/z):

  m/z 73   [Si(CH3)3]+                  a TMS group is present                                      needs >= MIN_73 (10 %)
  m/z 75   [(CH3)2Si=OH]+ (and 30Si)    backs m/z 73: other m/z 73 ions (C4H9O+, C3H5O2+) have no silicon   needs >= MIN_75 (1 %)
  m/z 147  [(CH3)2Si=O-Si(CH3)3]+       at least two TMS groups (suggested, see below); backs m/z 73   needs >= MIN_147 (2 %)
  m/z 115  [Si(CH3)2C(CH3)3]+           a silicon ion of TBDMS (with 73 and 75, the corroboration of M-57)   needs >= MIN_TBDMS_SI (2 %)
  M-15     loss of CH3 from a Si-CH3    confirms the molecular ion M (the usual strongest M ion)     needs >= MIN_M_ION (1 %)
  M-31     loss of OCH3                 MeOx (methoxime) evidence                                    needs >= MIN_M31 (2 %) and >= M31_OF_M15 x I(M-15)
  M-57     loss of C(CH3)3             TBDMS: the tert-butyl group goes; often base-peak-like       needs >= MIN_M57_TBDMS (20 %), > I(M-15), and a silicon ion
  M-90     loss of HOSi(CH3)3 (TMSOH)   shown for the user only; no decision reads it (see the measurements)
  M-180    loss of two TMSOH            shown for the user only; no decision reads it (see the measurements)

M is the molecular ion m/z of the DERIVATIVE the user typed; without it only m/z 73, 75, 115 and 147 are read. The thresholds are working values for
unit-mass EI spectra, not published cut-offs; they are constants here so one place changes them.

Inference (analyze(), used when the derivatization is Auto):
  TBDMS       M given, M-57 >= 20 % and stronger than M-15, and one of m/z 73, 75, 115 >= 2 %
  TMS         m/z 73 >= 10 % AND (m/z 75 >= 1 % or m/z 147 >= 2 %); or, when m/z 73 cannot be seen (the lowest peak of the spectrum is above m/z 73,
              scans that start at m/z 85 are common), m/z 147 >= MIN_147_BLIND (5 %). M-15, M-90 and M-180 are not evidence for TMS
  MeOx + TMS  TMS evident, M given and M-31 as above. TMS evident without M-31: MeOx stays "undecided" (the answer is TMS; the workflow tries both
              TMS and MeOx + TMS for compounds with a C=O and keeps the better match), because M-31 is weak or absent in most MeOx spectra
  None        no TMS and no TBDMS evidence
The lower bound of the TMS count is 1 with TMS evidence and 2 when m/z 147 is >= 2 %. The bound is a suggestion (see below); the exact count of a candidate
comes from its derivative structure (derivatization.derivatize), not from the spectrum.

How the numbers were chosen. Measured on the local MassBank 2026.03 EI index (13,477 spectra, 70 eV, mostly sector-instrument data from Japan), each
spectrum put through the workflow's own cleanup first. Labels: TMS ether = the entry's SMILES carries O/N/S-Si(CH3)3 (n = 441); name-only TMS = no SMILES,
"TMS" in the name (49); GC-TOF = a parent SMILES with a derivatized spectrum, from the GC-TOF metabolite sets of Kazusa, RIKEN, Osaka University and
GL Sciences, restricted to parents that have a silylation site (1,114; 76 % of them are recorded from above m/z 73); TBDMS = SMILES with Si(CH3)2C(CH3)3
(44); negative = SMILES without silicon, no derivative wording, highest peak <= M+30 (10,369; 116 spectra that are plainly derivatives under a parent
SMILES were taken out; "M known" = the compound's own nominal mass). Organosilicon compounds with only Si-C bonds (762) are neither class: about half of
them pass the TMS test, which no marker can change.

                                            before                              after (the constants above)
  TMS ether recognised                      82.5 % (no M), 88.4 % (M known)     81.0 % (M known or not)
  name-only TMS recognised                  69.4 %                              63.3 %
  GC-TOF recognised                         87.8 %                              82.9 %   (scans from above m/z 73: 84.4 % -> 78.1 %)
  negative called TMS (false positive)      6.8 % (no M), 11.1 % (M known)      2.7 %    (no SMILES: 4.1 % -> 3.2 %)
  negative called TMS, scans above m/z 73   11.6 % (n = 335)                    6.9 %    (3.3 % at 10 %, with 70 % of those GC-TOF recognised)
  TBDMS recognised (M known)                40.9 %                              40.9 %
  negative called TBDMS (M known)           7.5 %                               2.3 %
  MeOx + TMS called, TMS evident, with C=O  21.1 % (GC-TOF)                     14.2 %
  MeOx + TMS called, TMS evident, no C=O    8.7 % (GC-TOF), 10.6 % (TMS ether)  3.7 %, 4.1 %
Limits of these numbers: the library is mostly sector-instrument data; the GC-TOF labels are inferred from the parent structure; "negative" only means
"not labelled as a derivative"; only 44 TBDMS spectra exist, mostly synthetic intermediates (the real MTBSTFA metabolite spectra would show M-57 more
often). They are a first estimate for the rates, not a promise for one instrument.
Findings that changed the rules: (1) the old second route, m/z 147 with M-15, added 4.3 points of false positives (M-15 >= 1 % is in 33 % of
underivatized spectra and 86 % of TMS ethers) for 6 % more recognised ethers: dropped. (2) M-90 and M-180 do not separate the classes (>= 1 %: M-90 in
20 % of TMS ethers, 18 % of GC-TOF, 16 % of underivatized spectra; M-180 in 7 %, 9 %, 8 %), so they neither raise the count bound nor corroborate; they are
only shown. (3) m/z 147 >= 2 % means two or more TMS in 73 % of the TMS ethers and 95 % of the GC-TOF compounds, but mono-TMS compounds show it too (17 %
and 34 % of them), so the bound is a suggestion. (4) M-31 is a weak MeOx marker (14 % of the GC-TOF compounds with a C=O); it only tips the call when clear.
"""
from __future__ import annotations

TOLERANCE = 0.5                 # Da, unit-mass EI
MZ_TMS, MZ_TMS_OH, MZ_TMS2, MZ_TBDMS = 73, 75, 147, 115
MIN_73, MIN_75, MIN_147 = 10.0, 1.0, 2.0     # % of the base peak
MIN_147_BLIND = 5.0             # % of the base peak: m/z 147 alone, where m/z 73 is outside the recorded range
MIN_M_ION = 1.0                 # % of the base peak: M, M-15, M-90, M-180
MIN_M31, M31_OF_M15 = 2.0, 0.1  # % of the base peak, and the share of the M-15 intensity that M-31 must also reach
MIN_M57_TBDMS = 20.0            # % of the base peak
MIN_TBDMS_SI = 2.0              # % of the base peak: one of TBDMS_SI_IONS must reach it for a TBDMS call
TBDMS_SI_IONS = (MZ_TMS, MZ_TMS_OH, MZ_TBDMS)
TMS_NOMINAL, TBDMS_NOMINAL, MEOX_NOMINAL = 72, 114, 29       # nominal mass added by Si(CH3)3 (C3H8Si), Si(CH3)2C(CH3)3 (C6H14Si), C=O -> C=N-OCH3 (CH3N)
LOSSES = (("M-15", 15), ("M-31", 31), ("M-57", 57), ("M-90", 90), ("M-180", 180))

AUTO_NONE, AUTO_TMS, AUTO_MEOX_TMS, AUTO_TBDMS = "None", "TMS", "MeOx + TMS", "TBDMS"
MEOX_YES, MEOX_UNDECIDED, MEOX_NA = "yes", "undecided", "n/a"
_TMS_FAMILY = (AUTO_TMS, AUTO_MEOX_TMS)

MEANING = {"73": "TMS group", "75": "silicon ion, backs 73", "147": "two or more TMS groups (suggested)", "115": "TBDMS silicon ion", "M": "molecular ion",
           "M-15": "CH3 loss, confirms M", "M-31": "OCH3 loss, MeOx", "M-57": "tert-butyl loss, TBDMS",
           "M-90": "TMSOH loss, shown only", "M-180": "2 x TMSOH loss, shown only"}


def relative_intensity(peaks, mz, tolerance=TOLERANCE):
    """Largest intensity within +-tolerance of ``mz``, as % of the largest intensity of the spectrum (0.0 when there is no such peak)."""
    top = max((i for _, i in peaks), default=0.0)
    if not top > 0:
        return 0.0
    return max((i for m, i in peaks if abs(m - mz) <= tolerance), default=0.0) / top * 100.0


def _marker(marker_id, label, mz, rel, minimum, needs_m=False):
    """One row of the marker table: present is None while M is missing for an M-related ion."""
    present = None if rel is None else rel >= minimum
    return {"id": marker_id, "label": label, "mz": mz, "rel": rel, "min": minimum, "present": present, "needs_m": needs_m, "meaning": MEANING[marker_id]}


def analyze(peaks, molecular_ion=None) -> dict:
    """Read the marker ions of ``peaks`` [(mz, intensity)] (any scale) -> dict:

      markers        [{"id", "label", "mz", "rel" (% of base, None without M), "min", "present" (True | False | None without M), "needs_m", "meaning"}]
                     in this order: 73, 75, 147, 115, M, M-15, M-31, M-57, M-90, M-180
      blind_73       True when the lowest peak is above m/z 73, so its absence says nothing
      tms            "yes" | "no": TMS evidence (m/z 73 with 75 or 147; or m/z 147 >= MIN_147_BLIND where 73 cannot be seen; see the module docstring)
      tms_min        suggested lower bound of the TMS count (0 without TMS evidence)
      tbdms          True when M-57 >= MIN_M57_TBDMS, stronger than M-15, and a silicon ion (73, 75, 115) is >= MIN_TBDMS_SI
      meox           "yes" (M-31 present with TMS) | "undecided" (TMS, M-31 not seen or no M) | "n/a" (no TMS evidence)
      derivatization the inferred derivatization: "TMS" | "MeOx + TMS" | "TBDMS" | "None"
      m_check        "ok" (M or one of its losses seen) | "above spectrum" (not seen, and M-15 is above the highest peak) | "not seen" | None (no M)
      reasons        short sentences saying what the inference rests on
      warnings       sentences for the user (a molecular ion that does not fit the spectrum)
    """
    peaks = [(float(m), float(i)) for m, i in peaks if i > 0]
    m = None if molecular_ion is None else float(molecular_ion)

    def rel_at(loss):
        return None if m is None else relative_intensity(peaks, m - loss)

    rows = [_marker("73", f"m/z {MZ_TMS}", MZ_TMS, relative_intensity(peaks, MZ_TMS), MIN_73),
            _marker("75", f"m/z {MZ_TMS_OH}", MZ_TMS_OH, relative_intensity(peaks, MZ_TMS_OH), MIN_75),
            _marker("147", f"m/z {MZ_TMS2}", MZ_TMS2, relative_intensity(peaks, MZ_TMS2), MIN_147),
            _marker("115", f"m/z {MZ_TBDMS}", MZ_TBDMS, relative_intensity(peaks, MZ_TBDMS), MIN_TBDMS_SI),
            _marker("M", "M", m, None if m is None else relative_intensity(peaks, m), MIN_M_ION, True)]
    for name, loss in LOSSES:
        minimum = {"M-31": MIN_M31, "M-57": MIN_M57_TBDMS}.get(name, MIN_M_ION)
        rows.append(_marker(name, name, None if m is None else m - loss, rel_at(loss), minimum, True))
    by_id = {r["id"]: r for r in rows}
    seen = lambda key: bool(by_id[key]["present"])          # noqa: E731 - None (no M) counts as not seen
    rel = lambda key: by_id[key]["rel"]                     # noqa: E731

    low_mz = min((x for x, _ in peaks), default=0.0)
    blind_73 = low_mz > MZ_TMS + TOLERANCE                  # m/z 73 is outside the recorded range, so its absence says nothing
    si_ion = max((rel(str(k)) for k in TBDMS_SI_IONS), default=0.0)
    tbdms = seen("M-57") and rel("M-57") > rel("M-15") and si_ion >= MIN_TBDMS_SI
    backed_73 = seen("73") and (seen("75") or seen("147"))
    blind_147 = blind_73 and rel("147") >= MIN_147_BLIND
    tms = backed_73 or blind_147
    meox_31 = seen("M-31") and rel("M-31") >= M31_OF_M15 * rel("M-15")
    tms_min = (2 if seen("147") else 1) if tms else 0
    reasons, warnings = [], []
    if tbdms:
        derivatization, meox = AUTO_TBDMS, MEOX_NA
        reasons.append(f"M-57 is {rel('M-57'):.0f} % of the base peak and stronger than M-15 (tert-butyl loss), with a silicon ion (m/z 73, 75 or 115) at {si_ion:.0f} %")
    elif tms:
        meox = MEOX_YES if meox_31 else MEOX_UNDECIDED
        derivatization = AUTO_MEOX_TMS if meox == MEOX_YES else AUTO_TMS
        if backed_73:
            backing = [f"m/z 75 ({rel('75'):.0f} %)" if seen("75") else "", f"m/z 147 ({rel('147'):.0f} %)" if seen("147") else ""]
            reasons.append(f"m/z 73 is {rel('73'):.0f} % of the base peak, backed by " + " and ".join(x for x in backing if x))
        else:
            reasons.append(f"m/z 73 is outside the recorded range and m/z 147 is {rel('147'):.0f} % of the base peak")
        if meox == MEOX_YES:
            reasons.append("M-31 (OCH3 loss) is present")
        elif m is None:
            reasons.append("MeOx cannot be decided without the molecular ion")
        else:
            reasons.append("M-31 is not seen, so MeOx is undecided")
    else:
        derivatization, meox = AUTO_NONE, MEOX_NA
        if blind_73:
            first = (f"m/z 73 is outside the recorded range and m/z 147 is {rel('147'):.1f} % of the base peak (needs {MIN_147_BLIND:g} %)" if rel("147")
                     else "m/z 73 is outside the recorded range and m/z 147 is absent")
        elif seen("73"):
            first = f"m/z 73 is {rel('73'):.1f} % of the base peak but m/z 75 ({rel('75'):.1f} %) and m/z 147 ({rel('147'):.1f} %) do not back it"
        else:
            first = f"m/z 73 is {rel('73'):.1f} % of the base peak (below {MIN_73:g} %)"
        reasons.append(first + ("" if m is not None else "; give the molecular ion to look for TBDMS"))
    if tms and seen("147") and tms_min == 2:
        reasons.append("m/z 147 present (suggests two or more TMS groups)")

    m_check = None
    if m is not None:
        top_mz = max((x for x, _ in peaks), default=0.0)
        if any(seen(k) for k in ("M", "M-15", "M-31", "M-57", "M-90", "M-180")):
            m_check = "ok"
        elif m - 15 > top_mz + TOLERANCE:
            m_check = "above spectrum"
            warnings.append(f"The molecular ion m/z {m:g} is above the highest peak of the spectrum (m/z {top_mz:g}); check that value.")
        else:
            m_check = "not seen"
            warnings.append(f"None of M, M-15, M-31, M-57, M-90 or M-180 reaches its minimum for the molecular ion m/z {m:g}; "
                            "check that value (TMS derivatives usually show M-15).")
    return {"markers": rows, "blind_73": blind_73, "tms": "yes" if tms else "no", "tms_min": tms_min, "tbdms": tbdms, "meox": meox, "derivatization": derivatization,
            "m_check": m_check, "reasons": reasons, "warnings": warnings}


def check_selection(analysis, selected):
    """Warnings when the derivatization the user picked ("None" | "TMS" | "MeOx + TMS" | "TBDMS") is at odds with the markers. The choice is kept."""
    inferred, out = analysis["derivatization"], []
    by_id = {r["id"]: r for r in analysis["markers"]}
    if selected == AUTO_NONE and (analysis["tms"] == "yes" or analysis["tbdms"]):
        out.append(f"You chose None, but the markers point to {inferred} ({'; '.join(analysis['reasons'][:1])}). Your choice is kept.")
    elif selected in _TMS_FAMILY and analysis["tbdms"]:
        out.append(f"You chose {selected}, but M-57 points to TBDMS. Your choice is kept.")
    elif selected in _TMS_FAMILY and analysis["tms"] == "no" and not analysis["blind_73"]:
        if by_id["73"]["present"]:
            out.append(f"You chose {selected}, but m/z 73 ({by_id['73']['rel']:.1f} %) is not backed by m/z 75 or 147 (a TMS spectrum shows at least "
                       f"{MIN_75:g} % at m/z 75 or {MIN_147:g} % at m/z 147). Your choice is kept.")
        else:
            out.append(f"You chose {selected}, but m/z 73 is only {by_id['73']['rel']:.1f} % of the base peak (a TMS spectrum shows at least "
                       f"{MIN_73:g} %). Your choice is kept.")
    elif selected == AUTO_TBDMS and not analysis["tbdms"]:
        if by_id["M-57"]["present"] is False:
            out.append(f"You chose TBDMS, but M-57 is only {by_id['M-57']['rel']:.1f} % of the base peak (TBDMS shows at least {MIN_M57_TBDMS:g} %). Your choice is kept.")
        elif analysis["tms"] == "yes" and by_id["73"]["present"]:
            out.append("You chose TBDMS, but the markers look like TMS (m/z 73 is strong; give the molecular ion to test M-57). Your choice is kept.")
    return out


def state_of_spectrum(peaks, molecular_ion=None) -> str:
    """Derivatization state a reference spectrum itself shows: "TMS" | "TBDMS" | "None" (TMS and MeOx + TMS both read as TMS here)."""
    found = analyze(peaks, molecular_ion)["derivatization"]
    return AUTO_TMS if found == AUTO_MEOX_TMS else found


def family(mode):
    """The comparable family of a derivatization name: "TMS" (TMS, MeOx + TMS), "TBDMS", or "None"."""
    return AUTO_TMS if mode in _TMS_FAMILY else AUTO_TBDMS if mode in (AUTO_TBDMS, "MeOx + TBDMS") else AUTO_NONE
