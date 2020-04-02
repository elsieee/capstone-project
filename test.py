import bw
# x= [cycle gi0 go0 gi1 go1 gi2 go2]

def optimization():

	cost_list=[]
	for cycle in range(0,120):
		print("\n")
		print("cycle=",cycle)
		print("\n")
	
		for gi0 in range(0,cycle):
			for go0 in range(0,cycle):
				for gi1 in range(0,cycle):
					for go1 in range(0,cycle):
						for gi2 in range(0,cycle):
							for go2 in range(0,cycle):
								a = bw.Artery(cycle)
								a.addIntersection(bw.Intersection('i2', cycle, gi0, go0, 0.))
								a.addSegment(1000, 25, 35)
								a.addIntersection(bw.Intersection('i3', cycle, gi1, go1, 0.))
								a.addSegment(1000, 25, 35)
								a.addIntersection(bw.Intersection('i4', cycle, gi2, go2, 0.))

								bstar, bbarstar = a.optimize_pretimed_lp()
								print(bastar,bbarstar)
								print("\n")
								cost = bstar + bbstar
								cost_list.append(cost)

	max_value=max(cost_list)

	print(max_value)


	# for inter in a.intersection:
	# 	print(f"{inter.name}:{inter.absoffseto},{inter.absoffseti}")

	# print(bstar, bbarstar)

if __name__=='__main__':
	optimization()

#for b,c,d,e,f,g,h in zip(cycle gi0 go0 gi1 go1 gi2 go2): 
# ((cycle,gi0, go0, gi1, go1, gi2, go2) for cycle in [0,600] for gi0 in [0,cycle] for go0 in [0,cycle]
# 	for gi1 in [0,cycle] for go1 in [0,cycle] for gi2 in [0,cycle] for go2 in [0,cycle]):

######################################
#cycle = 600           # cycle in seconds