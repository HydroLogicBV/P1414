## Branches
class Branches:
    ruwheidhoog = 23
    ruwheidlaag = 23
    typeruwheid = "Bos en Bijkerk"


## Bridges
class Bridges:
    intreeverlies = 0.5
    lengte = 1
    ruwheid = 75.0
    typeruwheid = "StricklerKn"
    uittreeverlies = 0.7


## Culverts
class Culverts:
    breedteopening = 0.5
    gesloten = 1
    hoogtebinnenonderkantbene = -10
    hoogtebinnenonderkantbov = -10
    hoogteopening = 0.5
    intreeverlies = 0.6
    lengte = 1
    ruwheid = 75.0
    typeruwheid = "StricklerKn"
    uittreeverlies = 0.8
    vormkoker = 1


## Pumps
class Pumps:
    doelvariabele = "waterstand"
    maximalecapaciteit = 0
    peil_marge = 0.1


## Weirs
class Weirs:
    afvoercoefficient_stuw = 1
    afvoercoefficient_opening = 0.85
    hoogstedoorstroombreedte = None
    hoogstedoorstroomhoogte = 20
    laagstedoorstroombreedte = 0.1
    laagstedoorstroomhoogte = 20
    overlaatonderlaat = "Overlaat"
    soortregelbaarheid = 1
    soortstuw = 11
    vormopening = 3