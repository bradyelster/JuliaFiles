function time_to_mix_juice(juice)
    juice == "Pure Strawberry Joy" ? 0.5 :
    juice == "Green Garden" ? 1.5 :
    juice == "Energizer" ? 1.5 :
    juice == "Tropical Island" ? 3.0 :
    juice == "All or Nothing" ? 5.0 :
    2.5
end

function limes_to_cut(needed, limes)
    wedge_count = Dict("small" => 6, "medium" => 8, "large" => 10)
    wedges = 0
    limes_cut = 0

    for lime in limes
        wedges >= needed && break
        wedges += wedge_count[lime]
        limes_cut += 1
    end

    return limes_cut
end

function order_times(orders)
    [time_to_mix_juice(order) for order in orders]
end


function remaining_orders(time_left, orders)
    for (i, order) in enumerate(orders)
        time_needed = time_to_mix_juice(order)

        if time_left <= 0
            # No time left to START this order
            return orders[i:end]
        end

        # She starts this order and finishes it (even if it goes over time)
        time_left -= time_needed
    end

    # All orders completed
    return []
end