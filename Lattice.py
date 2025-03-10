"""
Author: tzs820, s240670, Emil Henningsen

Compute position and interrelative displacement vectors of different lattice constructions
"""

from utils import *

class Lattice:

    #Declarations:
    
    N = None
    d = None

    positions = None
    displacements = None
    polarizations = None

    def __init__(self, N : int = None, d : float = None, pos : ndarray = None, disp : ndarray = None, pola : ndarray = None):
        """
        TODO: Description
        """
        self.positions = pos
        self.displacements = disp
        self.polarizations = pola
        self.N = N
        self.d = d

    def clearLat(self) -> None:
        self.positions = None
        self.displacements = None
        self.polarizations = None
        self.N = None
        self.d = None

    """Get methods:"""

    def getN(self) -> int:
        return self.N
    
    def getd(self) -> float:
        return self.d

    def getPositions(self) -> ndarray:
        return self.positions
    
    def getDisplacements(self) -> ndarray:
        return self.displacements
    
    def getPolarizations(self) -> ndarray:
        return self.polarizations
    
    """Set methods:"""

    def setN(self, N : int) -> None:
        self.N = N
    
    def setd(self, d : float)-> None:
        self.d = d
        
    def setPositions(self, pos : ndarray)-> None:
        self.positions = pos

    def setDisplacements(self, disp : ndarray)-> None:
        self.displacements = disp
    
    def setPolarizations(self, pola : ndarray)-> None:
        self.polarizations = pola

    """Helper functions:"""

    def fillPola(self, N, dir : ndarray) -> ndarray:
        """
        TODO: Description
        """
        from numpy import zeros
        polarizations = zeros((N,3))
        for i in range(N):
            polarizations[i,:] = dir
        return polarizations

    def fillDisplacements(self, N : int, pos : ndarray) -> ndarray:
        """
        TODO: Description
        """
        from numpy import zeros

        rij = zeros((N, N, 3))                                       #Array of interrelative displacement vectors between lattice sites.
        for i in range(N):                                              #Fill the array
            for j in range(i, N):
                rij[i, j] = pos[i] - pos[j]
                rij[j, i] = -rij[i, j]                               #this array is symmetric by construction.

        return rij

    """Construct different standard lattices:"""

    def linlat(self, N : int, d : float, direction : ndarray = None, polarizations : ndarray = None, center : bool = True)-> None:
        """
        --- Parameters:

         - N: integer, number of lattice sites
         - d: float, distance between sites
         - direction (optional): array (3) of floats, representing lattice vector.
            if none given, defaults to x direction.
         - polarizations (optional): array (N x 3) of floats, representing each dipoles' polarization direction.
            if none given and internal container empty, defaults to z direction.
         - center (bool) - center the lattice around (0,0,0). 

        --- Return:

         - pos: array (N x 3) of floats, position vectors for each site
         - rij: array (N x N x 3) of floats, interrelative displacement vectors
        """
        from numpy import zeros, linspace, array, full_like, arange
        from utils import ex, ez

        #If already init, re-init.
        self.clearLat()

        #if no direction is supplied, set standard direction to x:
        if direction is None:
            direction = ex
        
        #if no polarization is supplied and not already set, set standard to z.
        if polarizations is None and self.polarizations is None:
            polarizations = self.fillPola(N, ez)                #Orient all dipoles along z
        elif polarizations.shape == (3,):
            dir = polarizations
            polarizations = self.fillPola(N, dir)
            #print("filled pola")

        pos = zeros((N, 3))                                          #Array of position vectors of the lattice sites (atoms).
        #ticks = linspace(-(N*d)/2, (N*d)/2, N)                       #displacement vectors in lattice LEGACY CALCULATION
        ticks = arange(0, N, 1) * d
        if center:
            ticks = ticks - (N-1)*d/2
        #ticks = []                                             #Brug denne udfyldning af ticks, for at sammenligne med Adrians resultater.
        #for i in range(N):
        #    ticks.append(i*d)
        for i in range(len(ticks)):
            pos[i, :] = direction * ticks[i]

        rij = self.fillDisplacements(N, pos)

        #Store in internal containers:
        self.setN(N)
        self.setd(d)
        self.setPositions(pos)
        self.setDisplacements(rij)
        self.setPolarizations(polarizations)

    def twopiece(self, N : int, d : float, theta : float, direction : ndarray = None, polarizations : ndarray = None) -> None:
        """
        Linear lattice broken in two pieces with angle theta in [0,2pi) on corner.
        TODO: Description
        """
        from utils import ex
        from numpy import array, cos, sin, concatenate, dot
        
        #if lattice is already initialized, clear out:
        self.clearLat()

        if direction is None:
            direction = ex

        first_dir = direction
        second_dir = array([cos(theta), 0, sin(theta)]) #break into z-direction by default. TODO: Hardcode

        #if N is even, break at half+1 lattice site
        if N % 2 == 0:
            N_first = int(N/2)
        else:
            N_first = int(N/2 + 1)
        N_second = N - N_first

        #Set polarizations:
        if polarizations is None:   #by default x-dir
            first_pola = ex
            second_pola = ex
        elif polarizations.shape == (N,3):
            first_pola = polarizations[0:N_first, :]
            second_pola = polarizations[N_first+1:-1, :]
        else:
            first_pola = polarizations
            second_pola = polarizations
        #else if pola is just vector, do nothing

        #Construct first linlat piece and store
        self.linlat(N_first, d, first_dir, first_pola)
        first_pos = self.getPositions()
        first_pola = self.getPolarizations()
        #Replace:
        for i in range(N_first):
            first_pos[i,:] -= N_first/2 * d * first_dir

        #Second:
        self.linlat(N_second, d, second_dir, second_pola)
        second_pos = self.getPositions()
        second_pola = self.getPolarizations()
        #Replace:
        for i in range(N_second):
            second_pos[i,:] += (N_second/2 + 1) * d * second_dir

        #Join lattices together:
        pos = concatenate([first_pos, second_pos], axis=0)
        disp = self.fillDisplacements(N, pos)
        pola = concatenate([first_pola, second_pola], axis=0)

        self.setPositions(pos)
        self.setDisplacements(disp)
        self.setPolarizations(pola)
        self.setN(N)
        self.setd(d)

    def circlelat(self, N : int, d : float, distance_measure : circle_measures = "inter", 
                  std_polarization : circle_polarizations = "other", polarizations : ndarray = None) -> None:
        """
        TODO: Description

        if std_polarization == "other", polarizations has to be supplied (either a single vec (3,) or full array (N,3)).
        """
        from numpy import zeros, pi, cos, sin, sqrt
        from utils import ez

        #If already init, clear lattice:
        self.clearLat()

        #do cirlce:
        def fillCircle(N : int, r : float, angle : float) -> ndarray:
            """
            TODO: Description
            """
            pos = zeros((N,3))
            for i in range(N):
                #fill positions for dipoles in xy-plane
                pos[i,0] = r * cos(i * angle)
                pos[i,1] = r * sin(i * angle)
            return pos
        b = 0                                           #inter-dipole distance
        r = 0                                           #radius
        angle = 2*pi / N
        match(distance_measure):
            case "inter":
                b = d
                r = b / sqrt(2 * (1 - cos(angle)))      #Cosine relation
            case "radius":
                r = d
                b = r * sqrt(2 * (1 - cos(angle)))      #Cosine relation
        pos = fillCircle(N, r, angle)

        #Displacements:
        rij = self.fillDisplacements(N, pos)

        #Do polarization:
        pola = None
        match(std_polarization):
            case "other":
                if polarizations is None and self.polarizations is None:
                    #Raise warning and choose standard direction
                    warn(self.warning_pola_none)
                    pola = self.fillPola(N, ez)
                elif polarizations.shape == (3,):
                    #provided polarization is a single vector, fill entire array.
                    pola = self.fillPola(N, polarizations)
                else:
                    #Last case, the entire array is supplied:
                    pola = polarizations
            case "inwards":
                pola = zeros((N,3))
                for i in range(N):
                    x = -cos(angle*i)
                    y = -sin(angle*i)
                    pola[i,0] = x
                    pola[i,1] = y
            case "radial alternating":
                pola = zeros((N,3))
                for i in range(N):
                    x = (-1)**i * cos(angle*i)
                    y = (-1)**i * sin(angle*i)
                    pola[i,0] = x
                    pola[i,1] = y
            case "outwards":
                pola = zeros((N,3))
                for i in range(N):
                    x = cos(angle*i)
                    y = sin(angle*i)
                    pola[i,0] = x
                    pola[i,1] = y
            case "azimuthal":
                pola = zeros((N,3))
                for i in range(N):
                    x = cos(angle*i + pi/2)
                    y = sin(angle*i + pi/2)
                    pola[i,0] = x
                    pola[i,1] = y
            case "azimuthal alternating":
                pola = zeros((N,3))
                for i in range(N):
                    x = (-1)**i * cos(angle*i + pi/2)
                    y = (-1)**i * sin(angle*i + pi/2)
                    pola[i,0] = x
                    pola[i,1] = y

        #Store in internal containers:
        self.setN(N)
        self.setd(d)
        self.setPositions(pos)
        self.setDisplacements(rij)
        self.setPolarizations(pola)