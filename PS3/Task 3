using Plots
# 1) #This is  a tuple that contians paraneters. we use par.x to use its contents if one uses just 'x' it wont work
par = (
    J  = 60,
    γ  = 2.0,
    γb = 1.0,
    β  = 0.96,
    r  = 0.04167,   
    θ  = 0.5,
    ā  = 2.0,       
    ȳ  = 1.0,       
    n_a = 500
)



# utility u(c) = c^(1-γ)/(1-γ)
function u(c, par)
    return c^(1 - par.γ) / (1 - par.γ)
end

#dla γb=1 -> log
function bequest(a, par)
    if par.γb == 1.0
        return par.θ * log(a + par.ā)
    else
        return par.θ * (a + par.ā)^(1 - par.γb) / (1 - par.γb)
    end
end


function Labor_income(par)
    y = zeros(par.J)
    for j in 1:par.J
        if j <= 40
            y[j] = par.ȳ * (0.8 + 0.02*j)
        else
            y[j] = par.ȳ * 0.3
        end
    end
    return y
end



function death_prob(par)
    pi_death = zeros(par.J)
    for i in 1:par.J
        pi_death[i] = min(0.0005*1.14^i, 1)
    end
    pi_death[par.J] = 1
    return pi_death
end

#if we consume nothing then 0+a'=(1+r)a+yj 
function find_a_max(par,y,a1=0.0)
    a_next=0.0
    a=a1
    for i in 1:par.J
        a_next=(1+par.r)*a+y[i]
        a=a_next
    end
    return a_next
end

#find the max value
y=Labor_income(par)
a_max=find_a_max(par,y)
par = merge(par, (a_max = a_max,)) #add to par


# create the grid
function get_a_grid(par)
    a_grid=ones(par.n_a)
    for i in 1:par.n_a
        a_grid[i]=par.a_max*((i-1)/(par.n_a-1))^2
    end
    return a_grid
end
a_grid=get_a_grid(par) #add to par
par=merge(par,(a_grid=a_grid,))

#============2===============#
#i write lots of comments for myslef as not to get lost this is hard actaully
function solve_model(par,y,pi_death)
    V = zeros(par.n_a, par.J)      # Value function
    pol_c = zeros(par.n_a, par.J)  # how much to consume given your assets
    pol_ap = zeros(par.n_a, par.J) # how much to save when agent is x years old
    #terminal period 
    for (i,a) in enumerate(par.a_grid) #for every possible asset
        cash = (1 + par.r) * a + y[par.J] #total money we have
        max_v=-Inf
        for (j, ap) in enumerate(par.a_grid) #loop trough ALL POSSIBLE SAVING CHOICES AVAILIBLE RN 
            if ap<=cash 
                c=cash-ap #cash= right side of budget constraint so ap=a'
                if c>0
                    value=u(c,par)+par.β*bequest(ap,par) #total utility
                    if value>max_v 
                        max_v=value
                        pol_c[i, par.J] = c
                        pol_ap[i, par.J] = ap
                    end
                end
            end
        end
        V[i,par.J]=max_v
    end
    
    for j in (par.J-1):-1:1 #loop backwards trough years
        for (i, a) in enumerate(par.a_grid)
            cash = (1 + par.r) * a + y[j]
            max_v = -Inf
            for (j_ap, ap) in enumerate(par.a_grid) #this is similar to what i've done before but now we have expected future value cause we dont imiedietly die
                if ap <= cash
                    c = cash - ap
                    if c > 0
                        continuation = (1 - pi_death[j]) * V[j_ap, j+1] + pi_death[j] * bequest(ap, par) #we dont die*value function at next year + we die and we leave  a bequest
                        value = u(c, par) + par.β * continuation #total value comsuption now + contuantion discounted 
                        if value > max_v
                            max_v = value
                            pol_c[i, j] = c              #this is the same as in terminal state
                            pol_ap[i, j] = ap            #
                        end
                    end
                end
            end
            V[i, j] = max_v
        end
    end
    
    return V, pol_c, pol_ap
end

pi_death = death_prob(par)
V, pol_c, pol_ap = solve_model(par, y, pi_death)

#PLOTTING PART 3 AND 4
ages = [20, 30, 40, 50, 60]

plot(par.a_grid, pol_c[:, 20], label = "Age 20")
plot!(par.a_grid, pol_c[:, 30], label = "Age 30")
plot!(par.a_grid, pol_c[:, 40], label = "Age 40")
plot!(par.a_grid, pol_c[:, 50], label = "Age 50")
plot!(par.a_grid, pol_c[:, 60], label = "Age 60")
title!("Consumption Policy")

plot(par.a_grid, pol_ap[:, 20], label = "Age 20")
plot!(par.a_grid, pol_ap[:, 30], label = "Age 30")
plot!(par.a_grid, pol_ap[:, 40], label = "Age 40")
plot!(par.a_grid, pol_ap[:, 50], label = "Age 50")
plot!(par.a_grid, pol_ap[:, 60], label = "Age 60")
title!("Savings Policy")

plot(par.a_grid, V[:, 20], label = "Age 20")
plot!(par.a_grid, V[:, 30], label = "Age 30")
plot!(par.a_grid, V[:, 40], label = "Age 40")
plot!(par.a_grid, V[:, 50], label = "Age 50")
plot!(par.a_grid, V[:, 60], label = "Age 60")
title!("Value Function")

consumption_path = zeros(60)
asset_path = zeros(60)
asset_path[1] = 0.0

for age in 1:60
    current_assets = asset_path[age]
    closest_index = argmin(abs.(par.a_grid .- current_assets))
    consumption_path[age] = pol_c[closest_index, age]
    if age < 60
        asset_path[age + 1] = pol_ap[closest_index, age]
    end
end

plot(1:60, consumption_path, label = "Consumption")
plot!(1:60, y, label = "Income")
plot!(1:60, asset_path, label = "Assets")
title!("Lifecycle")

#now i need to resolve, it shouldnt be hard cause only y_bar changes
#im going to brute force it
#change pars to pars_i literally re do the whole thing and plot
#this is not going to be pretty or efficient



#first one
par1 = merge(par, (ȳ = 0.5,))
y1 = Labor_income(par1)
a_max1 = find_a_max(par1, y1)
par1 = merge(par1, (a_max = a_max1,))
a_grid1 = get_a_grid(par1)
par1 = merge(par1, (a_grid = a_grid1,))
V1, pol_c1, pol_ap1 = solve_model(par1, y1, pi_death) #pi_death does not depend on ȳ it is ocnstant for all three versions 

consumption_path1 = zeros(60)
asset_path1 = zeros(60)
asset_path1[1] = 0.0
for age in 1:60
    current_assets = asset_path1[age]
    closest_index = argmin(abs.(par1.a_grid .- current_assets))
    consumption_path1[age] = pol_c1[closest_index, age]
    if age < 60
        asset_path1[age + 1] = pol_ap1[closest_index, age]
    end
end



#second one
par2 = merge(par, (ȳ = 1.0,))
y2 = Labor_income(par2)
a_max2 = find_a_max(par2, y2)
par2 = merge(par2, (a_max = a_max2,))
a_grid2 = get_a_grid(par2)
par2 = merge(par2, (a_grid = a_grid2,))
V2, pol_c2, pol_ap2 = solve_model(par2, y2, pi_death)

consumption_path2 = zeros(60)
asset_path2 = zeros(60)
asset_path2[1] = 0.0
for age in 1:60
    current_assets = asset_path2[age]
    closest_index = argmin(abs.(par2.a_grid .- current_assets))
    consumption_path2[age] = pol_c2[closest_index, age]
    if age < 60
        asset_path2[age + 1] = pol_ap2[closest_index, age]
    end
end




#tird one
par3 = merge(par, (ȳ = 2.0,))
y3 = Labor_income(par3)
a_max3 = find_a_max(par3, y3)
par3 = merge(par3, (a_max = a_max3,))
a_grid3 = get_a_grid(par3)
par3 = merge(par3, (a_grid = a_grid3,))
V3, pol_c3, pol_ap3 = solve_model(par3, y3, pi_death)

consumption_path3 = zeros(60)
asset_path3 = zeros(60)
asset_path3[1] = 0.0
for age in 1:60
    current_assets = asset_path3[age]
    closest_index = argmin(abs.(par3.a_grid .- current_assets))
    consumption_path3[age] = pol_c3[closest_index, age]
    if age < 60
        asset_path3[age + 1] = pol_ap3[closest_index, age]
    end
end



 #plot them

plot(1:60, consumption_path1, label = "(0.5)")
plot!(1:60, consumption_path2, label = "(1.0)")
plot!(1:60, consumption_path3, label = "(2.0)")
title!("Consumption")

plot(1:60, asset_path1, label = "(0.5)")
plot!(1:60, asset_path2, label = "(1.0)")
plot!(1:60, asset_path3, label = "(2.0)")
title!("Assets")

savings_rate1 = (asset_path1[2:60] - asset_path1[1:59]) ./ y1[1:59]
savings_rate2 = (asset_path2[2:60] - asset_path2[1:59]) ./ y2[1:59]
savings_rate3 = (asset_path3[2:60] - asset_path3[1:59]) ./ y3[1:59]

plot(1:59, savings_rate1, label = "0.5)")
plot!(1:59, savings_rate2, label = "(1.0)")
plot!(1:59, savings_rate3, label = "(2.0)")
title!("Savings Rate")

#BEWARE!
#btw i have a very big mess in vs code and some things broke and uplouding them via vs code to github is broken forsome reaosn idk why but im too tired to wrestle with this 
#i will do it trough github online so i will literally delet all previous code and upload new one ;-D

#NOTE FOR TASK 6 literlaly do the same as above but do params5, params6 etc(tough there will be lots of them cause 3 dofferent ȳ and two θ) for a lazy but working solution tough it nt be pretty to read
#NEED to discuss teh result of increaisng y/y_bar
#possible bug fixes needed its to late rn to check
#resolvng the modle can be (should be) made better 
