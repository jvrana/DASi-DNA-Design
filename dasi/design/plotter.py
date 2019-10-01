import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# TODO: plot primers
class Plotter:
    def __init__(self, length, starts, ends, types, names, costs, title=""):
        segments = zip(starts, ends, types, names, costs)
        self.segments = list(segments)
        self.length = length

        self.yspacer = 0.3
        self.border = 1
        self.text_spacer = 0.2
        self.title = title
        self.segment_styles = {
            "fragment": {"color": "black", "linewidth": 3},
            "synthesis": {"color": "red", "linewidth": 3},
            "overhang": {"color": "blue", "linewidth": 1},
        }

    def plot(self):

        Y = len([f for f in self.segments if f[2] != "overhang"])

        fig = plt.figure(figsize=(5, self.yspacer * Y))

        ax = fig.gca()
        ax.set_xlim(-self.border, self.length + self.border)
        ax.set_ylim(-self.border, Y + self.border)

        # set border and axis
        ax.axes.get_xaxis().set_visible(True)
        ax.axes.get_yaxis().set_visible(False)
        ax.set_xlabel("bp")
        ax.set_frame_on(False)
        ax.set_title(self.title)
        y = Y
        for i, seg in enumerate(self.segments):
            x1, x2, segtype, name, cost = seg
            if "overlap" in name or "overlap" in segtype:
                segtype = "overhang"
            elif "synthesis" in name or "synthesis" in segtype:
                segtype = "synthesis"
            else:
                segtype = "fragment"
            cost = "{:.1f}".format(float(cost))
            name = "{} - {}".format(Y - y, name)
            segstyle = self.segment_styles[segtype]
            if x1 < x2:
                line = Line2D([seg[0], seg[1]], [y, y], zorder=1, **segstyle)

                # add cpst
                xtext = (seg[0] + seg[1]) / 2.0
                ax.text(
                    xtext,
                    y + self.text_spacer,
                    str(cost),
                    horizontalalignment="center",
                    verticalalignment="bottom",
                )

                # add label
                ax.text(
                    self.length,
                    y,
                    name,
                    horizontalalignment="left",
                    verticalalignment="center",
                )
                ax.add_line(line)
                y -= 1
            else:
                if segtype == "overhang":
                    if i == len(self.segments) - 1:
                        y2 = Y
                    else:
                        y2 = y + 1
                    line1 = Line2D([x1, x2], [y, y2], zorder=0, **segstyle)
                    line2 = Line2D([x2, x1], [y, y2], zorder=0, **segstyle)

                    xtext = (x1 + x2) / 2.0
                    ax.text(
                        xtext,
                        y - 2 * self.text_spacer,
                        str(cost),
                        horizontalalignment="center",
                        verticalalignment="top",
                    )

                    ax.add_line(line1)
                    ax.add_line(line2)
                else:
                    line1 = Line2D([x1, self.length], [y, y], zorder=1, **segstyle)
                    line2 = Line2D([0, x2], [y, y], zorder=1, **segstyle)

                    # add cost
                    if self.length - x1 > x2:
                        xtext = (self.length - x1) / 2.0
                    else:
                        xtext = x2 / 2.0
                    ax.text(
                        xtext,
                        y + self.text_spacer,
                        str(cost),
                        horizontalalignment="center",
                        verticalalignment="bottom",
                    )

                    # add label
                    ax.text(
                        self.length,
                        y,
                        name,
                        horizontalalignment="left",
                        verticalalignment="center",
                    )

                    ax.add_line(line1)
                    ax.add_line(line2)
        plt.show()
