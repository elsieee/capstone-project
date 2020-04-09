def optimization_3inter():
    
    max_cost = 0
    for cycle in range(60,121,30):
        print("\n")
        print("cycle=", cycle)
        print("\n")
        
        for gi0 in range(0,cycle, int(cycle/3)):
            for go0 in range(0,cycle, int(cycle/3)):
                for gi1 in range(0,cycle, int(cycle/3)):
                    for go1 in range(0,cycle, int(cycle/3)):
                        for gi2 in range(0,cycle, int(cycle/3)):
                            for go2 in range(0,cycle, int(cycle/3)):
                                
                                
                                a = bw.Artery(cycle)
                                a.addIntersection(bw.Intersection('i2', cycle, gi0, go0, 0))
                                a.addSegment(200, 25, 35)
                                a.addIntersection(bw.Intersection('i3', cycle, gi1, go1, 0))
                                a.addSegment(200, 25, 35)
                                a.addIntersection(bw.Intersection('i4', cycle, gi2, go2, 0))

                                bstar, bbarstar = a.optimize_pretimed_lp()
                                cost = bstar + bbarstar
                                if cost >= max_cost:
                                    max_cost = cost                                    
                                    optimized_parameters = [cycle, gi0, go0, gi1, go1, gi2, go2]
                                    

    print("\n")
    print("max_cost =", max_cost)
    print("\n")
    print (optimized_parameters)
    return max_cost


if __name__ == '__main__':
    optimization_3inter()
    
