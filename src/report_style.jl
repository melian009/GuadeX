"""
    Shared style for the report figures.

Every report figure is inserted in the LaTeX reports at a fixed physical width
(about 16 cm), so what determines the apparent text size is the ratio

    font size / image logical width.

A figure that is 3000 logical pixels wide with 20 pt text renders its labels at
roughly 3 pt once scaled to 16 cm, which is unreadable.  This file therefore
enforces two conventions:

* **Large fonts relative to the canvas.** The report figures deliberately use
  large font sizes (base 30, axis labels 30, tick labels 26, legend 26, panel
  letters 34, annotations 26) and compact figure canvases, so the text occupies
  a much larger fraction of each image than in the earlier versions.
* **No descriptive figure titles.** The explanation lives in the LaTeX
  captions.  Multi-panel figures are identified with bold letters ``a``, ``b``,
  ``c`` ... placed in the axis title slot, left-aligned, so the caption can
  refer to individual panels.

`save_report_figure` also trims the uniform white border that Makie leaves
around a layout; trimming is what makes the remaining canvas map onto the
printed width, which further increases the apparent font size.
"""

const REPORT_BASE_FONTSIZE = 30
const REPORT_AXIS_LABEL_FONTSIZE = 30
const REPORT_TICK_FONTSIZE = 26
const REPORT_LEGEND_FONTSIZE = 26
const REPORT_PANEL_LETTER_FONTSIZE = 34
const REPORT_ANNOTATION_FONTSIZE = 26

"""
    report_theme!(; base=REPORT_BASE_FONTSIZE)

Apply the shared report-figure theme.  Returns `nothing`; the theme is global
state, so call it at the top of every report-figure function.
"""
function report_theme!(; base::Real=REPORT_BASE_FONTSIZE)
    set_theme!(Theme(
        fontsize = base,
        Axis = (
            xlabelsize = REPORT_AXIS_LABEL_FONTSIZE,
            ylabelsize = REPORT_AXIS_LABEL_FONTSIZE,
            xticklabelsize = REPORT_TICK_FONTSIZE,
            yticklabelsize = REPORT_TICK_FONTSIZE,
            titlesize = REPORT_PANEL_LETTER_FONTSIZE,
            titlefont = :bold,
            titlealign = :left,
        ),
        Legend = (
            labelsize = REPORT_LEGEND_FONTSIZE,
            titlefontsize = REPORT_LEGEND_FONTSIZE,
        ),
        Colorbar = (
            labelsize = REPORT_AXIS_LABEL_FONTSIZE,
            ticklabelsize = REPORT_TICK_FONTSIZE,
        ),
    ))
    return nothing
end

"""
    panel_letter!(ax, letter)

Write a bold panel letter (for example `"a"` or `"(a)"`) above the left edge of
`ax`, in the title slot so it never overlaps the plotted data.
"""
function panel_letter!(ax, letter::AbstractString)
    ax.title = letter
    ax.titlesize = REPORT_PANEL_LETTER_FONTSIZE
    ax.titlefont = :bold
    ax.titlealign = :left
    return ax
end

const _ANNOTATION_POSITIONS = Dict(
    :lt => (0.03, 0.97, (:left, :top)),
    :lb => (0.03, 0.04, (:left, :bottom)),
    :rt => (0.97, 0.97, (:right, :top)),
    :rb => (0.97, 0.04, (:right, :bottom)),
)

"""
    report_annotation!(ax, text; position=:lt)

Add a plain in-panel annotation (not a title) using the report annotation font
size.  `position` is one of `:lt`, `:lb`, `:rt`, `:rb`.
"""
function report_annotation!(ax, text::AbstractString; position::Symbol=:lt, color=:black)
    x, y, align = get(_ANNOTATION_POSITIONS, position, _ANNOTATION_POSITIONS[:lt])
    return text!(ax, x, y; space = :relative, text = text, align = align,
        fontsize = REPORT_ANNOTATION_FONTSIZE, color = color)
end

"""
    equal_panel_columns!(fig, weights...)

Force the layout columns of `fig` to the given relative widths.  Makie sizes
`Auto` columns from their content, which makes a column that holds a legend or a
colourbar different from a column of axes; giving every column an explicit
weight keeps the data panels the same size.
"""
function equal_panel_columns!(fig, weights::Real...)
    for (c, w) in enumerate(weights)
        colsize!(fig.layout, c, Auto(w))
    end
    return fig
end

"""
    trim_figure!(path)

Remove the uniform white border around a saved figure when ImageMagick is
available.  A trimmed image maps its content onto the printed column width, so
the text appears larger.  No-op (and never an error) when `magick` is missing.
"""
function trim_figure!(path::AbstractString)
    magick = Sys.which("magick")
    magick === nothing && return path
    isfile(path) || return path
    tmp = path * ".trim.png"
    try
        run(`$magick $path -trim +repage $tmp`)
        isfile(tmp) && mv(tmp, path; force=true)
    catch err
        # Never fail a figure build because of the optional trim step.
        isfile(tmp) && rm(tmp; force=true)
        @warn "could not trim figure; keeping untrimmed image" path exception=err
    end
    return path
end
