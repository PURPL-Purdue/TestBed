if __name__ == "__main__":
    ThrustStart = 12
    ThrustMax = 500
    Pcstart = 100
    PcMax = 300
    OFstart = .5
    OFEnd = 4
    FireTime = 2

    while Pcstart <= PcMax:
        while ThrustStart <= ThrustMax:
            while OFstart <= OFEnd:
                CEA(ThrustStart, OFstart, Pcstart)

                OFflowChecker(OFstart, mdot, FireTime)




    



def OFflowChecker(OF, mdot, FireTime):
    {

    
    }

def plotter(OF):
    {


    }

def CEA(thrust, OF, Pc):
    {
        

    }