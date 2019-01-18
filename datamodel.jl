using DataFrames, CSV, Ipopt, JuMP

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
  for i in nodes
    if i == length(nodes)
      push!(arcs, (i,1))
    else
      push!(arcs, (i,i+1))
    end
  end
  data[:arcs] = arcs

  return data

end

data = create_data_model()

VRP = Model(solver=IpoptSolver())

#########################################################
#   Mechanical Power consumption of vehicle_capacity    #
#########################################################
velocity = 70*rand()   #Choose a random speed between 0 and 70
mass_vehicle = 500     #Mass of vehicle
mass_person = 80*rand() #Choose the weight of one person
gravity = 9.8          #acceleration due to gravity
duration = 3*rand()    #time to travel between arc i and i+1
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
#nodes = range(1, data[ :num_locations])
@variable(VRP, x[i in nodes], lowerbound = 0, upperbound = 1, start = 1)
@variable(VRP, dur[i in nodes], lowerbound = 1)         #arrival time at node i
@variable(VRP, P_arc[i in nodes]) #power consumption after traversing from node i to i+1
@variable(VRP, F_arc[i in nodes]) #fuel rate over an arc
@variable(VRP, F_con[i in nodes]) #fuel consumption over arc i to i+1
@variable(VRP, acc[i in nodes])   #acceleration of vehicle from node i to i+1
@variable(VRP, ppload[i in nodes], lowerbound = 1) #number of people in each node
@variable(VRP, TotalCost, lowerbound=FV)

##################################
#           Constraints          #
##################################
#Power consumption on each arc
@constraint(VRP, PC[i in nodes],
        (P_arc[i] == velocity * (Fa + coeff2*cos(alpha) + mass_vehicle*gravity*sin(alpha) + mass_person*gravity*ppload[i])))
@constraint(VRP, FR[i in nodes],
        (F_arc[i] == (ξ/(k*ψ))*(k*N*D) + P_arc[i]/(η*ηtf))) #fuel rate for each arc
@constraint(VRP, FC[i in nodes],
        (F_con[i] == dur[i] * F_arc[i]))    #fuel consumption on each arc
@constraint(VRP, NoV,
        (sum(x[i] for i in nodes) ≤ data[:num_vehicles]))  #vehicle on each route is less than or equal to number of vehicles in fleet
@constraint(VRP, MaxCapacity,
        (sum(ppload[i] for i in nodes) ≤ data[:vehicle_capacity])) #restrict number of people to maximum capacity of vehicle
#@constraint(VRP, MaxTimePerNode,
#        ((dur[i] ≤ data[:wait][1][i] for i in nodes)))
        #doesn't seem to work yet  and I think it's better than the next constraint Maximum waiting time at each node
@constraint(VRP, MaxTime,
        (sum(dur[i] for i in nodes) ≤ sum(data[:wait][1][i] for i in nodes))) #Total time of traversing all nodes

@constraint(VRP, ObjectiveFunction,
 (TotalCost == sum(ζ*x[i]*F_con[i] for i in nodes)))  #Total cost of running x_vehicles over all nodes

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



# function prold(a)
#   b = 1
#   for i in range(1,a)
#     b = b*(i+1)
#   end
#   return b/(a+1)
# end
# prold(20)
