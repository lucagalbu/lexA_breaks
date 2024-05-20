package me.lucagalbusera.simplebreaks

case class Rates(
    mrnaCreation: Double,
    mrnaDecay: Double,
    breakRepair: Double,
    proteinFolding: Double,
    lexaUnbinding: Double,
    proteinCreation: Double
)

object Rates:
  // create a rates with the mrna creation set to have a specific slope in the peak distribution
  def createFromSlope(
      mrnaDecay: Double,
      breakRepair: Double,
      proteinFolding: Double,
      slope: Double = RecaSlope,
      lexaUnbinding: Double = LexaUnbinding,
      proteinCreation: Double = ProteinCreation
  ): Rates =
    val x = mrnaDecay / proteinFolding
    val alpha =
      if x == 1 then proteinCreation / math.E
      else proteinCreation * math.pow(x, -x / (x - 1))
    val exponent = math.exp(slope * alpha / CellVolume)
    val mrnaCreation = exponent / (1 - exponent) * breakRepair / CellVolume

    Rates(
      mrnaCreation = mrnaCreation,
      mrnaDecay = mrnaDecay,
      proteinFolding = proteinFolding,
      breakRepair = breakRepair,
      proteinCreation = proteinCreation,
      lexaUnbinding = lexaUnbinding
    )
