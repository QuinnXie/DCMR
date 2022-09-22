using PlotlyJS

function plotRunningTime(runningTime, start, alg::Collocation, MRMethod)
    nStart = length(start)
    nMRMethod = length(MRMethod)

    if (typeof(alg) <: Trapezoidal{Continuous})
        plot_title = "Time for Solving Dynamic Optimization Probelm (TRPuC)"
    elseif (typeof(alg) <: Trapezoidal{Discontinuous})
        plot_title = "Time for Solving Dynamic Optimization Probelm (TRPuD)"
    elseif (typeof(alg) <: HermiteSimpson{Continuous})
        plot_title = "Time for Solving Dynamic Optimization Probelm (HSSuC)"
    elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
        plot_title = "Time for Solving Dynamic Optimization Probelm (HSSuD)"
    else
        @error "Unexpected method and iput rule '$(alg)'"
    end

    pltRunTime = PlotlyJS.plot(
        [
            PlotlyJS.bar(name=start[1], x=MRMethod, y=runningTime[1:nStart:end-nStart]),
            PlotlyJS.bar(name=start[2], x=MRMethod, y=runningTime[2:nStart:end-nStart]),
            PlotlyJS.bar(name=start[3], x=MRMethod, y=runningTime[3:nStart:end-nStart]),
            PlotlyJS.bar(name=start[4], x=MRMethod, y=runningTime[4:nStart:end-nStart]),
            PlotlyJS.bar(name=start[1], x=[MRMethod[end]], y=runningTime[end-nStart+1:nStart:end], yaxis="y2"),
            PlotlyJS.bar(name=start[2], x=[MRMethod[end]], y=runningTime[end-nStart+2:nStart:end], yaxis="y2"),
            PlotlyJS.bar(name=start[3], x=[MRMethod[end]], y=runningTime[end-nStart+3:nStart:end], yaxis="y2"),
            PlotlyJS.bar(name=start[4], x=[MRMethod[end]], y=runningTime[end-nStart+4:nStart:end], yaxis="y2")
        ],
        Layout(
            xaxis_title_text="Mesh Refinement Method",
            yaxis_title_text="Running Time (s)",
            yaxis2=attr(
                title="Running Time (s) (only for MRBisectionOne)",
                overlaying="y",
                side="right"
            ),
            title=attr(
                text=plot_title,
                y=0.95,
                x=0.5,
                xanchor= "center",
                yanchor= "top"
            ),
            legend=attr(
                x=1,
                y=-0.3,
                yanchor="bottom",
                xanchor="right",
                orientation="h",
                bordercolor="Black",
                borderwidth=2
            ),
        )
    )
    display(pltRunTime)

    # touch("Cart_Pole_v4_plot\\runningTime\\warmStart.png")
    # io  = open("warmStart.png", "w")
    # PlotlyJS.savefig(
    #     io,
    #     # MIME"image/eps"(),
    #     pltRunTime,
    #     format="png"
    # )

    # PlotlyJS.savefig(pltRunTime, format="png")
    return pltRunTime
end

function plotRunningTime(runningTime, start, collocationMethod, meshRefi::MeshRefi)
    nStart = length(start)
    ncollocationMethod = length(collocationMethod)

    if typeof(meshRefi) <: MRMinimax
        MRType = "MRMinimax"
    elseif typeof(meshRefi) <: MREquidistributed
        MRType = "MREquidistributed"
    elseif typeof(meshRefi) <: MRBisectionMulti
        MRType = "MRBisectionMulti"
    elseif typeof(meshRefi) <: MRBisectionAll
        MRType = "MRBisectionAll"
    elseif typeof(meshRefi) <: MRBisectionOne
        MRType = "MRBisectionOne"
    else
        MRType = "$(typeof(meshRefi))"
    end

    plot_title = "Time for Solving Dynamic Optimization Probelm ($MRType)"

    pltRunTime = PlotlyJS.plot(
        [
            PlotlyJS.bar(name=start[1], x=collocationMethod, y=runningTime[1:nStart:end]),
            PlotlyJS.bar(name=start[2], x=collocationMethod, y=runningTime[2:nStart:end]),
            PlotlyJS.bar(name=start[3], x=collocationMethod, y=runningTime[3:nStart:end]),
            PlotlyJS.bar(name=start[4], x=collocationMethod, y=runningTime[4:nStart:end]),
            # PlotlyJS.bar(name=start[1], x=[MRMethod[end]], y=runningTime[end-nStart+1:nStart:end], yaxis="y2"),
            # PlotlyJS.bar(name=start[2], x=[MRMethod[end]], y=runningTime[end-nStart+2:nStart:end], yaxis="y2"),
            # PlotlyJS.bar(name=start[3], x=[MRMethod[end]], y=runningTime[end-nStart+3:nStart:end], yaxis="y2"),
            # PlotlyJS.bar(name=start[4], x=[MRMethod[end]], y=runningTime[end-nStart+4:nStart:end], yaxis="y2")
        ],
        Layout(
            xaxis_title_text="Mesh Refinement Method",
            yaxis_title_text="Running Time (s)",
            # yaxis2=attr(
            #     title="Running Time (s) (only for MRBisectionOne)",
            #     overlaying="y",
            #     side="right"
            # ),
            title=attr(
                text=plot_title,
                y=0.95,
                x=0.5,
                xanchor= "center",
                yanchor= "top"
            ),
            legend=attr(
                x=1,
                y=-0.3,
                yanchor="bottom",
                xanchor="right",
                orientation="h",
                bordercolor="Black",
                borderwidth=2
            ),
        )
    )
    display(pltRunTime)
    # PlotlyJS.savefig(pltRunTime, format="png")
    return pltRunTime
end

# alg = Collocation(Ipopt.Optimizer(), 4, Continuous(), "Trapezoidal")
# runningTime = medianTime[1:20]
# plotRunningTime(runningTime, start, alg::Collocation, MRMethod)

# alg = Collocation(Ipopt.Optimizer(), 4, Discontinuous(), "Trapezoidal")
# runningTime = medianTime[21:40]
# plotRunningTime(runningTime, start, alg::Collocation, MRMethod)

# alg = Collocation(Ipopt.Optimizer(), 2, Continuous(), "Hermite-Simpson")
# runningTime = medianTime[41:60]
# plotRunningTime(runningTime, start, alg::Collocation, MRMethod)

# alg = Collocation(Ipopt.Optimizer(), 2, Discontinuous(), "Hermite-Simpson")
# runningTime = medianTime[61:80]
# pltRunTime = plotRunningTime(runningTime, start, alg::Collocation, MRMethod)
# # PlotlyJS.savefig(pltRunTime, format="pdf")
# # savefig(pltRunTime, format="png")
# # plotRunningTimeNoBisectionOne(runningTime, start, alg, MRMethod)

# function plotRunningTimeNoBisectionOne(runningTime, start, alg::Collocation, MRMethod)
#     nStart = length(start)
#     nMRMethod = length(MRMethod)

#     if (typeof(alg) <: Trapezoidal{Continuous})
#         plot_title = "Trapezoidal Collocation with Continuous Input"
#     elseif (typeof(alg) <: Trapezoidal{Discontinuous})
#         plot_title = "Trapezoidal Collocation with Discontinuous Input"
#     elseif (typeof(alg) <: HermiteSimpson{Continuous})
#         plot_title = "Hermite-Simpson Collocation with Continuous Input"
#     elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
#         plot_title = "Hermite-Simpson Collocation with Discontinuous Input"
#     else
#         @error "Unexpected method and iput rule '$(alg)'"
#     end

#     pltRunTime = PlotlyJS.plot(
#         [
#             PlotlyJS.bar(name=start[1], x=MRMethod, y=runningTime[1:nStart:end]),
#             PlotlyJS.bar(name=start[2], x=MRMethod, y=runningTime[2:nStart:end]),
#             PlotlyJS.bar(name=start[3], x=MRMethod, y=runningTime[3:nStart:end]),
#             PlotlyJS.bar(name=start[3], x=MRMethod, y=runningTime[3:nStart:end]),
#             PlotlyJS.bar(name=start[4], x=MRMethod, y=runningTime[4:nStart:end]),
#             # PlotlyJS.bar(name=start[1], x=[MRMethod[end]], y=runningTime[end-nStart+1:nStart:end], yaxis="y2"),
#             # PlotlyJS.bar(name=start[2], x=[MRMethod[end]], y=runningTime[end-nStart+2:nStart:end], yaxis="y2"),
#             # PlotlyJS.bar(name=start[3], x=[MRMethod[end]], y=runningTime[end-nStart+3:nStart:end], yaxis="y2"),
#             # PlotlyJS.bar(name=start[4], x=[MRMethod[end]], y=runningTime[end-nStart+4:nStart:end], yaxis="y2")
#         ],
#         Layout(
#             xaxis_title_text="Mesh Refinement Method",
#             yaxis_title_text="Running Time (s)",
#             # yaxis2=attr(
#             #     title="Running Time (s) (only for MRBisectionOne)",
#             #     overlaying="y",
#             #     side="right"
#             # ),
#             title=attr(
#                 text=plot_title,
#                 y=0.95,
#                 x=0.5,
#                 xanchor= "center",
#                 yanchor= "top"
#             ),
#             legend=attr(
#                 x=1,
#                 y=-0.2,
#                 yanchor="bottom",
#                 xanchor="right",
#                 orientation="h",
#                 bordercolor="Black",
#                 borderwidth=2
#             ),
#         )
#     )
#     display(pltRunTime)
#     # PlotlyJS.savefig(pltRunTime)
# end
    # p1 = Array{GenericTrace{Dict{Symbol, Any}},1}(undef, nStart)
    # for i in 1:nStart
    #     p1[i] = PlotlyJS.bar(name=start[i], x=MRMethod, y=runningTime[i:nStart:end-4])
    # end

    # p2 = Array{GenericTrace{Dict{Symbol, Any}},1}(undef, nStart)
    # for i in 1:nStart
    #     p2[i] = PlotlyJS.bar(name=start[i], x=[MRMethod[end]], y=runningTime[end-nStart+i:nStart:end])
    # end
    # # p1 = PlotlyJS.bar(name=start[1], x=MRMethod, y=runningTime[1:nStart:end])
    # # p2 = PlotlyJS.bar(name=start[2], x=MRMethod, y=runningTime[2:nStart:end])
    # # p3 = PlotlyJS.bar(name=start[3], x=MRMethod, y=runningTime[3:nStart:end])
    # # p4 = PlotlyJS.bar(name=start[4], x=MRMethod, y=runningTime[4:nStart:end])

    # layout = Layout(
    #     # xaxis_title="Mesh Refinement Method",
    #     yaxis_title="Running Time (s)",
    #     title=attr(
    #         # text="Trapezoidal Collocation with Continuous Input",
    #         font_family="Times New Roman",
    #         font_color="black",
    #         y=1.1,
    #         x=0.5,
    #         xanchor= "center",
    #         yanchor= "top"
    #     ),
    #     font=attr(
    #         family="Times New Roman",
    #         size=17,
    #         color="black"
    #     ),
    #     legend_title_text="Initialization",
    #     legend=attr(
    #         # x=0,
    #         # y=1,
    #         title_font_family="Times New Roman",
    #         font=attr(
    #             family="Times New Roman",
    #             size=15,
    #             color="black"
    #         ),
    #         x=1,
    #         y=0.9,
    #         yanchor="bottom",
    #         xanchor="right",
    #         orientation="h",
    #         # bgcolor="LightSteelBlue",
    #         bordercolor="Black",
    #         borderwidth=2
    #     ),
    # )

    # if (typeof(alg) <: Trapezoidal{Continuous})
    #     layout.title[:text] = "Trapezoidal Collocation with Continuous Input"
    #     plot_title = "Trapezoidal Collocation with Continuous Input"
    # elseif (typeof(alg) <: Trapezoidal{Discontinuous})
    #     layout.title[:text] = "Trapezoidal Collocation with Discontinuous Input"
    #     plot_title = "Trapezoidal Collocation with Discontinuous Input"
    # elseif (typeof(alg) <: HermiteSimpson{Continuous})
    #     layout.title[:text] = "Hermite-Simpson Collocation with Continuous Input"
    #     plot_title = "Hermite-Simpson Collocation with Continuous Input"
    # elseif (typeof(alg) <: HermiteSimpson{Discontinuous})
    #     layout.title[:text] = "Hermite-Simpson Collocation with Discontinuous Input"
    #     plot_title = "Hermite-Simpson Collocation with Discontinuous Input"
    # else
    #     @error "Unexpected method and iput rule '$(alg)'"
    # end

    # plt1 = PlotlyJS.plot([p1[i] for i in 1:nStart], layout)

    # layout.legend[:orientation] = "v"
    # layout.legend[:y] = "0.73"
    # layout.yaxis_title = " "
    # plt2 = PlotlyJS.plot([p2[i] for i in 1:nStart], layout)
    # plt3 = [plt1 plt2]
    # relayout!(plt3, height=500, width=700, title_text=plot_title)
    # display(plt3)
