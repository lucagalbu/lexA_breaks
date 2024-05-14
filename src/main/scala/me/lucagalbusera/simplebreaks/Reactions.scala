package me.lucagalbusera.simplebreaks

import me.lucagalbusera.simplebreaks.RandomGenerator
import me.lucagalbusera.simplebreaks.RandomGenerator.randomGenerator

enum EnumReactions:
  case mrnaCreation, mrnaDecay, proteinCreation, proteinFolded,
    breakRepaired, lexaBinding, lexaUnbinding, None

type ListReactions = List[(EnumReactions, Double)]

case class NextReaction(
    reaction: EnumReactions,
    deltaTime: Double
)

object NextReaction:
  def computeNext(
      rates: Rates,
      cell: Cell,
      randomGenerator: RandomGeneratorTrait
  ) =
    val reactions = reactionsFromRates(rates = rates, cell = cell)
    val probabilities = normalize(reactions = reactions)
    NextReaction(
      reaction = findNextReaction(
        probabilities = probabilities,
        randomGenerator = randomGenerator
      ),
      deltaTime = findNextDeltaTime(
        reactions = reactions,
        randomGenerator = randomGenerator
      )
    )

  private def findNextReaction(
      probabilities: ListReactions,
      randomGenerator: RandomGeneratorTrait
  ): EnumReactions =
    val cumProbabilities = computeCumSum(reactions = probabilities)
    val rnd = randomGenerator.getUnifRandom()
    cumProbabilities.filter((reaction, value) => value >= rnd).head(0)

  private def findNextDeltaTime(
      reactions: ListReactions,
      randomGenerator: RandomGeneratorTrait
  ): Double =
    val sumRates = reactions.map((_, rate) => rate).sum
    randomGenerator.getExpRandom(rate = sumRates)

    /** Given a set of reaction rate constants, it compute the reaction rates by
      * adding the current cell state. E.g., if the mRNA production rate
      * constant is rm and the cell is repressed, then the mRNA production rate
      * will be 0. If the cell is not repressed, then the mRNA production rate
      * will be rm.
      *
      * @param rates
      *   Constant rates of the reactions.
      * @param cell
      *   Current state of the cell.
      * @return
      */
  private def reactionsFromRates(rates: Rates, cell: Cell): ListReactions =
    val mrnaCreation = if cell.repressed then 0 else rates.mrnaCreation
    val mrnaDecay = rates.mrnaDecay * cell.mrna
    val proteinFolding = rates.proteinFolding * cell.protein
    val breakRepair = if cell.break then rates.breakRepair else 0
    val proteinCreation = rates.proteinCreation * cell.mrna
    val lexaUnbinding = if cell.repressed then rates.lexaUnbinding else 0
    val lexaBinding = if cell.repressed then 0 else cell.lexa

    List(
      EnumReactions.mrnaCreation -> mrnaCreation,
      EnumReactions.mrnaDecay -> mrnaDecay,
      EnumReactions.breakRepaired -> breakRepair,
      EnumReactions.lexaBinding -> lexaBinding,
      EnumReactions.lexaUnbinding -> lexaUnbinding,
      EnumReactions.proteinCreation -> proteinCreation,
      EnumReactions.proteinFolded -> proteinFolding
    )

  private def normalize(reactions: ListReactions): ListReactions =
    val totalSum = reactions.map((_, value) => value).sum
    reactions.map((name, value) => name -> value / totalSum)

  private def computeCumSum(reactions: ListReactions): ListReactions =
    reactions.scan[(EnumReactions, Double)]((EnumReactions.None, 0))(
      (prev, cur) => (cur(0), prev(1) + cur(1))
    )
