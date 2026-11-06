import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
D = "results/figures/loci/"
def build(names, out):
    ims = [mpimg.imread(D + n + ".png") for n in names]
    h = [i.shape[0] / i.shape[1] for i in ims]
    fig, ax = plt.subplots(len(ims), 1, figsize=(7.2, 7.2 * sum(h)), gridspec_kw=dict(height_ratios=h, hspace=0.02))
    for k, (a, im) in enumerate(zip(ax, ims)):
        a.imshow(im); a.axis("off")
        a.text(-0.01, 1.0, "abc"[k], transform=a.transAxes, fontsize=14, fontweight="bold", va="top", ha="right")
    fig.savefig(out + ".png", dpi=300, bbox_inches="tight"); fig.savefig(out + ".pdf", dpi=300, bbox_inches="tight")
build(["short_RAI1_chr17_17752484", "product_ZMIZ1_chr10_79150279", "short_HM13_chr20_31547393"], D + "Fig6")
build(["nested_CHMP1B", "internal_ACIN1_chr14_23071564"], D + "FigS10")
