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
      slope: Double,
      mrnaDecay: Double,
      breakRepair: Double,
      proteinFolding: Double,
      lexaUnbinding: Double = LexaUnbinding,
      proteinCreation: Double = ProteinCreation
  ): Rates =
    val V = 2.57
    val x = mrnaDecay / proteinFolding
    val alpha =
      if x == 1 then proteinCreation / math.E
      else proteinCreation * math.pow(x, -x / (x - 1))
    val exponent = math.exp(slope * alpha / V)
    val mrnaCreation = exponent / (1 - exponent) * breakRepair / V

    Rates(
      mrnaCreation = mrnaCreation,
      mrnaDecay = mrnaDecay,
      proteinFolding = proteinFolding,
      breakRepair = breakRepair,
      proteinCreation = proteinCreation,
      lexaUnbinding = lexaUnbinding
    )
