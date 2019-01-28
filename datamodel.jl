using DataFrames, CSV, JuMP, Gurobi, Cbc, Ipopt

function create_data_model()
  #Creates the data for the example.
  data = Dict{}()
  # Array of distances between locations.
  distances = readtable("distances.csv")

  data[:distances] = distances
  data[:num_locations] = length(distances)
  data[:num_vehicles] = 3
  data[:vehicle_capacity] = 20
  data[:depot] = 1  #Take-off point

  #possible demands per node
  demands = readtable("demands.csv")
  data[:demands] = demands
  #maximum waiting time per node
  wait = readtable("max_waiting_time.csv")
  data[:wait] = wait

  #######################
  #      Parameters     #
  #######################
  #All nodes
  global nodes = range(1, data[ :num_locations])

  #all arcs/roads
  arcs = []
  global N_arcs = 0
  for i in nodes
    if i == length(nodes)
      break #push!(arcs, (i,1))
      #N_arcs += 1
    else
      push!(arcs, (i,i+1))
      N_arcs += 1
    end
  end
  data[:arcs] = arcs
  data[:N_arcs] = N_arcs

  return data

end

data = create_data_model()

#I changed from Ipopt because Ipopt has no support for integer or binary variables
VRP = Model(solver=GurobiSolver())

#########################################################
#   Mechanical Power consumption of vehicle_capacity    #
#########################################################
#avg_distance = sum(data[:distances][1])/length(data[:distances][1])
velocity = rand(30:55)   #Choose a random speed between 0 and 70
mass_vehicle = 500     #Mass of vehicle
mass_person = rand(45:100) #Choose the weight of one person
gravity = 9.8          #acceleration due to gravity
duration = rand(0.5:2.5)    #time to travel between arc i and i+1
acc = velocity/duration     #acceleration of vehicle
alpha = 0.08727        #gradient angle of the road/arc
coeff1 = 0.05          #Rolling friction coefficient for cars on moving on snow
coeff2 = 0.37          #aerodynamic drag coefficient for Ford Transit Custom Mk8
air_Density = 1.3      #Air density for Moscow winter climate
area_frontal = 120     #Frontal area of the vehicle
ζ = 12                 #cost of diesel fuel per km
ξ = 14.7               #fuel_air_mass_ratio
k = 45.5               #heating value of diesel fuel in MJ/Kg
N = 750                #Engine speed in RPM
D = 1.57               #engine displacement if the engine bore is 10 cm and the stroke is 5 cm with four cylinders,
ψ = 0.85               #afactor converting the fuel rate from grams per second to liters per second
η = 0.90               #efficiency parameter for diesel engines
ηtf = 0.85             #drive train efficiency.

#Power is given by product of force and velocity
Fa = 0.5*air_Density*area_frontal*coeff2*(velocity^2) #aerodynamic Resistance
Fg = mass_vehicle*acc*sin(alpha)                      #Gravitational force
Fr = coeff1*mass_vehicle*gravity*cos(alpha)           #Rolling resistance
FV = mass_vehicle*acc                                 #weight of vehicle

#Total Mechanical power required by one vehicle is
PW = (Fa+Fg+Fr+FV)*velocity

##################################
#          Variables             #
##################################
#Major decision variable which is 1 if vehivle i goes between arcs i & i+1
@variable(VRP, x[i in 1:N_arcs], start = 1, Bin)
@variable(VRP, dur[i in 1:N_arcs], lowerbound = 1)         #arrival time at node i
@variable(VRP, P_arc[i in 1:N_arcs]) #power consumption after traversing from node i to i+1
@variable(VRP, F_arc[i in 1:N_arcs]) #fuel rate over an arc
@variable(VRP, F_con[i in 1:N_arcs]) #fuel consumption over arc i to i+1
@variable(VRP, acc[i in 1:N_arcs])   #acceleration of vehicle from node i to i+1
@variable(VRP, ppload[i in nodes], lowerbound = 0, Int) #number of people carried in each node
@variable(VRP, wait_dur[i in nodes])
@variable(VRP, TotalCost, lowerbound=FV)

##################################
#           Constraints
##################################
#Power consumption on each arc
@constraint(VRP, PowerConsumption[i in 1:N_arcs],
        (P_arc[i] == velocity * (Fa + coeff2*cos(alpha) + mass_vehicle*gravity*sin(alpha) + mass_person*gravity*ppload[i])))
@constraint(VRP, FuelRate[i in 1:N_arcs],
        (F_arc[i] == (ξ/(k*ψ))*(k*N*D) + P_arc[i]/(η*ηtf))) #fuel rate for each arc
#The following is an equality quadratic constraint which isn't supported by Gurobi
#But when I tried converting to inequality, I got another error about
#"Q" matrix not being positive definite
@constraint(VRP, FuelConsumption[i in 1:N_arcs],
        (F_con[i] == dur[i] * F_arc[i]))    #fuel consumption on each arc
@constraint(VRP, duration_arc[i in 1:N_arcs],
        (dur[i] ≤ data[:distances][1][i+1]/velocity))
@constraint(VRP, No_of_Vehicles,
        (sum(x[i] for i in N_arcs) ≤ data[:num_vehicles]))  #vehicle on each route is less than or equal to number of vehicles in fleet
@constraint(VRP, MaxCapacity,
        (sum(ppload[i] for i in nodes) ≤ data[:vehicle_capacity])) #restrict number of people to maximum capacity of vehicle#@constraint(VRP, MaxTimePerNode,
@constraint(VRP, NodeCapacity[i in nodes],
        (sum(ppload[i] for i in nodes) ≤ data[:vehicle_capacity]))
@constraint(VRP, Waittime[i in 1:N_arcs],
        (wait_dur[i] ≤ data[:wait][1][i])) #wait time at each node is lower than the assigned wait time per node
@constraint(VRP, MaxTime,
        (sum(dur[i] for i in N_arcs) ≤ sum(data[:wait][1][i] for i in nodes))) #Total time of traversing all nodes

#Same issue here with quadratic equality constraint
#I thought about setting all "x"s to 1?
@constraint(VRP, ObjectiveFunction,
 (TotalCost ≥ sum(ζ*x[i]*F_con[i] for i in 1:N_arcs)))  #Total cost of running x_vehicles over all nodes

##################################
#           Objective            #
##################################
@objective(VRP, Min, TotalCost)

print(VRP)

solve(VRP)

println("Objective Value: ", getobjectivevalue(VRP))
println("x: ", getvalue(x))
println("Duration: ", getvalue(dur))
println("Passengers: ", getvalue(ppload))
println("TotalCost: ", getvalue(TotalCost))
