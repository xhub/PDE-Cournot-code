# SPDX-License-Identifier: MIT


function consumer_surplus_discrete(qout::Vector{Float64}, times::Vector{Float64}, econ_data)

	middle_vals = collect(zip(times, qout))[2:end-1]
	Δt = times[2] - times[1]

	return Δt * sum((get_price(econ_data, tk, 0.)         - get_price(econ_data, tk, qk)) * qk/2 for (tk, qk) in middle_vals)
	      + Δt/2 * ((get_price(econ_data, times[1], 0.)   - get_price(econ_data, times[1], qout[1])) * qout[1]/2
					     +  (get_price(econ_data, times[end], 0.) - get_price(econ_data, times[end], qout[end])) * qout[end]/2)
end

function producer_surplus(qin::Vector{Float64}, qout::Vector{Float64}, times::Vector{Float64}, econ_data)

	middle_vals = collect(zip(times, qin, qout))[2:end-1]
	Δt = times[2] - times[1]

#	println("DEBUG: $(collect(get_price(econ_data, tk, qoutk) for (tk,qink,qoutk) in middle_vals))")
#println("DEBUG: $(collect(qoutk for (tk,qink,qoutk) in middle_vals))")
	return Δt * sum(get_price(econ_data, tk, qoutk)*qoutk - get_cost(econ_data, tk, qink) * qink for (tk,qink,qoutk) in middle_vals)
			+ Δt/2 * (get_price(econ_data, times[1], qout[1])*qout[1] - get_cost(econ_data, times[1], qin[1]) * qin[1]
						 +  get_price(econ_data, times[end], qout[end])*qout[end] - get_cost(econ_data, times[end], qin[end]) * qin[end])

end
