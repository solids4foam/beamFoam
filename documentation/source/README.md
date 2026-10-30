# Sources of the PDFs in documentation/

- `monolithicCouplingPlan.html` → `../monolithicCouplingPlan.pdf`
- `outerCorrectorsAndNewmark.html` → `../outerCorrectorsAndNewmark.pdf`, with figures
  `growth.png` and `newton.png` from `figs.py`

`model.py` is the linear body–spring model of the loop and monolithic updates
(amplification matrices). `order.py` checks its order of accuracy. Run them with
any Python that has numpy and matplotlib.

To render a PDF (macOS, from this directory):

    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" --headless=new \
        --disable-gpu --no-pdf-header-footer \
        --print-to-pdf=../outerCorrectorsAndNewmark.pdf "file://$PWD/outerCorrectorsAndNewmark.html"
