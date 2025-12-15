"""
    plot_catchment_network(g::AbstractGraph; labels=nothing, coordinates=nothing)

Visualizes a river subcatchment network.
- `g`: The directed graph object.
- `labels`: (Optional) Vector of names for the sites.
- `coordinates`: (Optional) Vector of Point2f or Tuple (x,y) for fixed geographic positions.
"""
function plot_catchment_network(g::AbstractGraph;
  labels=nothing,
  coordinates=nothing)

  # 1. Set up the Figure
  f = Figure(size=(800, 600))
  ax = Axis(f[1, 1], title="River Subcatchment Network")

  # Hide axis decorations (grid, ticks, borders) for a cleaner map look
  hidedecorations!(ax)
  hidespines!(ax)

  # 2. Determine Layout
  # If coordinates are provided, use them.
  # Otherwise, use a tree-friendly layout (BuchheimWalker) if available, or Spring.
  if !isnothing(coordinates)
    layout_func = coordinates
  else
    # BuchananWalker is excellent for river trees (roots at top/side)
    layout_func = NetworkLayout.spring(g)
  end

  # 3. Define Aesthetics
  # Color nodes based on degree (hubs are darker) or a fixed water blue
  node_color_map = :azure2
  edge_color_map = :slategray

  # 4. Plot the Graph
  p = graphplot!(ax, g;
    layout=layout_func,

    # Node styling
    node_size=20,
    node_color=node_color_map,
    node_strokewidth=2,
    node_strokecolor=:blue,

    # Edge styling (Arrows are critical for stream flow)
    arrow_show=true,
    arrow_size=15,
    edge_width=1.5,
    edge_color=edge_color_map,

    # Labels (if provided)
    nlabels=labels,
    nlabels_align=(:center, :bottom),
    nlabels_textsize=11,
    nlabels_distance=15
  )

  # 5. Return the figure object
  return f
end