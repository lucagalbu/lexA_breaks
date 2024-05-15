import org.scalamock.scalatest.MockFactory
import org.scalatest.funsuite.AnyFunSuite
import me.lucagalbusera.simplebreaks.NextReaction
import me.lucagalbusera.simplebreaks.ListReactions
import me.lucagalbusera.simplebreaks.Rates
import me.lucagalbusera.simplebreaks.Cell
import me.lucagalbusera.simplebreaks.RandomGeneratorTrait
import org.scalamock.matchers.Matchers
import me.lucagalbusera.simplebreaks.EnumReactions

class NextReactionTests extends AnyFunSuite with MockFactory with Matchers {
  test("The correct reaction is found when only two reactions can happen") {
    // With the followin test cell and reaction constants,
    // only breakRepair or lexaUnbinding can happen. The probabilities are 0.2 and 0.8
    val testCell = Cell(
      time = 0,
      lexa = 200,
      break = true,
      repressed = true,
      mrna = 0,
      protein = 0
    )

    val testRates = Rates(
      mrnaCreation = 1,
      mrnaDecay = 2,
      breakRepair = 2,
      proteinFolding = 4,
      lexaUnbinding = 8,
      proteinCreation = 7
    )

    val mockRandomGenerator = mock[RandomGeneratorTrait]
    (mockRandomGenerator.getExpRandom).expects(10).returns(0.1)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.6)

    val nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )

    assert(
      nextReaction == NextReaction(
        reaction = EnumReactions.lexaUnbinding,
        deltaTime = 0.1
      )
    )
  }

  test(
    "The last reaction is correctly found"
  ) {
    // All reactions happen except lexa unbinding
    val testCell = Cell(
      time = 0,
      lexa = 10,
      break = true,
      repressed = false,
      mrna = 1,
      protein = 1
    )

    val testRates = Rates(
      mrnaCreation = 15,
      mrnaDecay = 25,
      breakRepair = 5,
      proteinFolding = 15,
      lexaUnbinding = 10, // lexa binding is given by the lexa in the cell
      proteinCreation = 30
    )

    val mockRandomGenerator = mock[RandomGeneratorTrait]
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.99)

    val nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )

    assert(
      nextReaction == NextReaction(
        reaction = EnumReactions.proteinFolded,
        deltaTime = 0.1
      )
    )
  }

  test(
    "The first reaction is correctly found"
  ) {
    // All reactions happen except lexa unbinding
    val testCell = Cell(
      time = 0,
      lexa = 10,
      break = true,
      repressed = false,
      mrna = 1,
      protein = 1
    )

    val testRates = Rates(
      mrnaCreation = 15,
      mrnaDecay = 25,
      breakRepair = 5,
      proteinFolding = 15,
      lexaUnbinding = 10, // lexa binding is given by the lexa in the cell
      proteinCreation = 30
    )

    val mockRandomGenerator = mock[RandomGeneratorTrait]
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.00001)

    val nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )

    assert(
      nextReaction == NextReaction(
        reaction = EnumReactions.mrnaCreation,
        deltaTime = 0.1
      )
    )
  }

  test(
    "The correct reaction is found when a reaction has rate 0"
  ) {
    // Lexa unbinding doesn't happen. In this case the cumulative sums of probs have two
    // terms, lexaUnbinding and the next one, with the same cumultaive sum.
    val testCell = Cell(
      time = 0,
      lexa = 10,
      break = true,
      repressed = false,
      mrna = 1,
      protein = 1
    )

    val testRates = Rates(
      mrnaCreation = 15,
      mrnaDecay = 25,
      breakRepair = 10,
      proteinFolding = 10,
      lexaUnbinding = 0,
      proteinCreation = 30
    )

    val mockRandomGenerator = mock[RandomGeneratorTrait]
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.6)

    val nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )

    assert(
      nextReaction == NextReaction(
        reaction = EnumReactions.lexaBinding,
        deltaTime = 0.1
      )
    )
  }

  test(
    "Subsequent rates are correctly found, when the cell doesn't update"
  ) {
    val testCell = Cell(
      time = 0,
      lexa = 10,
      break = true,
      repressed = false,
      mrna = 1,
      protein = 1
    )

    val testRates = Rates(
      mrnaCreation = 15,
      mrnaDecay = 25,
      breakRepair = 10,
      proteinFolding = 10,
      lexaUnbinding = 0,
      proteinCreation = 30
    )

    val mockRandomGenerator = mock[RandomGeneratorTrait]
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)
    (mockRandomGenerator.getExpRandom).expects(100).returns(0.1)

    /* Rates cumulative probabilities are
    mrnaCreation -> 0.15, mrnaDecay -> 0.40, breakRepaired -> 0.50, lexaBinding -> 0.60,
    lexaUnbinding -> 0.60, proteinCreation -> 0.90, proteinFolded -> 1
     */
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.14)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.3)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.45)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.55)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.70)
    (() => mockRandomGenerator.getUnifRandom()).expects().returns(0.95)

    var nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )
    assert(nextReaction.reaction == EnumReactions.mrnaCreation)

    nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )
    assert(nextReaction.reaction == EnumReactions.mrnaDecay)

    nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )
    assert(nextReaction.reaction == EnumReactions.breakRepaired)

    nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )
    assert(nextReaction.reaction == EnumReactions.lexaBinding)

    nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )
    assert(nextReaction.reaction == EnumReactions.proteinCreation)

    nextReaction = NextReaction.computeNext(
      cell = testCell,
      rates = testRates,
      randomGenerator = mockRandomGenerator
    )
    assert(nextReaction.reaction == EnumReactions.proteinFolded)

  }
}
