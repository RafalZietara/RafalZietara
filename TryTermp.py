import sys
# sys.path.append('C:/Apps/SPPTS for Python/bin')

import shellthermo

solver = shellthermo.EquilibriumSolverBuilder().withDatabase('C:\Apps\SPPTS for Python\Pepper\experimental.thermodb').withComponents(["oH2", "pH2"]).withModel("CPA/SMIRK").withPhases("vll").build()

T = 298.15
stream = solver.isothermal(T, 1.0e5, [0.75, 0.25])
print(stream.T(), stream.P(), stream.numberOfPhasesPresent())
print(stream.T(), stream.P(), stream.density())


for i in range(stream.numberOfPhasesPossible()-1):
  if stream.isPhasePresent(i):
    p = stream.phase(i)
    print("Phase %d found, type %s" % (i, p.state()))
    print("Fraction: %s" % p.fraction())
    print("Composition:", p.componentFractions())
    print("Density:", p.density())
    print("Therm cond:", p.thermalConductivity())


try:
  # this should fail because there is no present phase at index 4
  print(stream.phase(4).fraction())
except Exception as err:
  print("Expected error occurred: %s" % err)