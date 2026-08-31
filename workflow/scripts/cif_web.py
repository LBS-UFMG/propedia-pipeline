"""Make an extracted per-pair mmCIF friendly to the website's 3D viewer.

Biopython `MMCIFIO` writes `_atom_site.label_atom_id` / `_atom_site.label_comp_id`
but NOT the `auth_` variants. 3Dmol.js (the viewer the site uses) reads the atom
name and residue name ONLY from `_atom_site.auth_atom_id` / `_atom_site.auth_comp_id`
— with no fallback to the `label_` columns — so without them the atom names come out
empty and 3Dmol cannot compute secondary structure (the cartoon renders as flat coil
instead of helices/sheets).

`add_auth_atom_comp` appends `auth_atom_id`/`auth_comp_id` columns to the `_atom_site`
loop, copying the `label_` values. That copy is exact and lossless: the atom NAME and
residue NAME are identical between the `auth` and `label` schemes — only `asym_id`
(chain) and `seq_id` (residue number) can ever differ, and those are untouched here.

The function is idempotent (a cif that already has `auth_atom_id` is returned
unchanged) and format-preserving (it only appends two fields per ATOM/HETATM row and
two header lines; every other line is passed through verbatim).
"""


def add_auth_atom_comp(text: str) -> str:
    if "_atom_site.auth_atom_id" in text:
        return text
    lines = text.split("\n")
    hdr_idx = [k for k, l in enumerate(lines) if l.strip().startswith("_atom_site.")]
    if not hdr_idx:
        return text                      # no atom_site loop -> nothing to do
    cols = [lines[k].strip() for k in hdr_idx]
    try:
        i_atom = cols.index("_atom_site.label_atom_id")
        i_comp = cols.index("_atom_site.label_comp_id")
    except ValueError:
        return text                      # unexpected layout -> leave untouched
    last_hdr = hdr_idx[-1]

    out = []
    for k, l in enumerate(lines):
        out.append(l)
        if k == last_hdr:
            out.append("_atom_site.auth_atom_id")
            out.append("_atom_site.auth_comp_id")

    res = []
    for l in out:
        t = l.strip()
        if t.startswith("ATOM ") or t.startswith("HETATM"):
            f = t.split()
            if len(f) > max(i_atom, i_comp):
                res.append(l + " " + f[i_atom] + " " + f[i_comp])
                continue
        res.append(l)
    return "\n".join(res)
