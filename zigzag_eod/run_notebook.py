"""Execute eod_zigzag_persistence.ipynb in-process and write outputs back.

nbconvert --execute needs to bind a socket for the ZMQ kernel, which this
sandbox forbids, so cells are exec'd in a shared namespace instead. Stdout,
errors, and figures are captured per cell and stored as notebook outputs, which
is what nbconvert would have produced.
"""
import io
import sys
import traceback
from contextlib import redirect_stdout

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nbformat as nbf

NB = "eod_zigzag_persistence.ipynb"
nb = nbf.read(NB, as_version=4)
ns = {"__name__": "__main__"}
failed = []

for i, cell in enumerate(nb.cells):
    if cell.cell_type != "code":
        continue
    buf = io.StringIO()
    cell.outputs = []
    try:
        with redirect_stdout(buf):
            exec(cell.source, ns)
    except Exception:
        tb = traceback.format_exc()
        cell.outputs.append(nbf.v4.new_output(
            "error", ename="Error", evalue="see traceback",
            traceback=tb.splitlines()))
        failed.append((i, tb))
        print(f"CELL {i} FAILED\n{tb}", file=sys.stderr)
        break
    finally:
        text = buf.getvalue()
        if text:
            cell.outputs.insert(0, nbf.v4.new_output("stream", name="stdout", text=text))
        # A cell that created a figure gets it attached as a display_data output.
        for num in plt.get_fignums():
            f = plt.figure(num)
            png = io.BytesIO()
            f.savefig(png, format="png", dpi=110, bbox_inches="tight")
            cell.outputs.append(nbf.v4.new_output(
                "display_data",
                data={"image/png": __import__("base64").b64encode(png.getvalue()).decode()},
                metadata={}))
            plt.close(f)
    cell.execution_count = i + 1
    print(f"cell {i:2d} ok", flush=True)

nbf.write(nb, NB)
print("FAILED" if failed else "ALL CELLS OK")
sys.exit(1 if failed else 0)
