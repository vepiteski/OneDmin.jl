# bmark_solvers  avec des OrderedDict

function my_bmark_solvers(solvers :: OrderedDict{Symbol,<: Any}, args...;
                       kwargs...)
  stats = OrderedDict{Symbol, DataFrame}()
  for (name, solver) in solvers
    @debug "running" name solver
    stats[name] = solve_problems(solver, args...; kwargs...)
  end
  return stats
end



# performance_profile  avec un OrderedDict
function my_performance_profile(stats::OrderedDict{Symbol,DataFrame}, cost::Function, args...; kwargs...)
  solvers = keys(stats)
  dfs = (stats[s] for s in solvers)
  P = hcat([cost(df) for df in dfs]...)
  performance_profile(P, string.(solvers), args...; kwargs...)
end




# profile_solvers  avec un OrderedDict
function my_profile_solvers(stats::OrderedDict{Symbol,DataFrame},
                            costs::Vector{<:Function},
                            costnames::Vector{String};
                            width::Int=400,
                            height::Int=400
                            )
  solvers = collect(keys(stats))
  dfs = (stats[solver] for solver in solvers)
  Ps = [hcat([cost(df) for df in dfs]...) for cost in costs]

  nprobs = size(stats[first(solvers)], 1)
  nsolvers = length(solvers)
  ncosts = length(costs)
  npairs = div(nsolvers * (nsolvers - 1), 2)
  colors = get_color_palette(:auto, nsolvers)

  # profiles with all solvers
  ps = [performance_profile(Ps[1], string.(solvers), palette=colors, title=costnames[1], legend=:bottomright)]
  nsolvers > 2 && xlabel!(ps[1], "")
  for k = 2 : ncosts
    p = performance_profile(Ps[k], string.(solvers), palette=colors, title=costnames[k], legend=false)
    nsolvers > 2 && xlabel!(p, "")
    ylabel!(p, "")
    push!(ps, p)
  end

  if nsolvers > 2
    ipairs = 0
    # combinations of solvers 2 by 2
    for i = 2 : nsolvers
      for j = 1 : i-1
        ipairs += 1
        pair = [solvers[i], solvers[j]]
        dfs = (stats[solver] for solver in pair)
        Ps = [hcat([cost(df) for df in dfs]...) for cost in costs]

        clrs = [colors[i], colors[j]]
        p = performance_profile(Ps[1], string.(pair), palette=clrs, legend=:bottomright)
        ipairs < npairs && xlabel!(p, "")
        push!(ps, p)
        for k = 2 : length(Ps)
          p = performance_profile(Ps[k], string.(pair), palette=clrs, legend=false)
          ipairs < npairs && xlabel!(p, "")
          ylabel!(p, "")
          push!(ps, p)
        end
      end
    end
    p = plot(ps..., layout=(1 + ipairs, ncosts), size=(ncosts * width, (1 + ipairs) * height))
  else
    p = plot(ps..., layout=(1, ncosts), size=(ncosts * width, height))
  end
  p
end

