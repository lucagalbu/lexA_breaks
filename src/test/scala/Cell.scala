import org.scalatest.funsuite.AnyFunSuite
import me.lucagalbusera.simplebreaks.Cell
import org.scalatest.funspec._
import org.scalatest.matchers.should._
import me.lucagalbusera.simplebreaks.NextReaction
import me.lucagalbusera.simplebreaks.EnumReactions
import me.lucagalbusera.simplebreaks.LexaDecay
import me.lucagalbusera.simplebreaks.LexaSteadyState

class CellTest extends AnyFunSuite with Matchers {
  test("LexA is updated correctly during a break") {
    var cell = Cell(time = 0, lexa = LexaSteadyState, break = true)

    for deltaTime <- Seq(0, 1, 1.5, 3, 3.5) do
      cell = Cell.updateCell(
        reaction =
          NextReaction(reaction = EnumReactions.None, deltaTime = deltaTime),
        cell = cell
      )
      val expectedLexa = LexaSteadyState * math.exp(-LexaDecay * cell.time)
      cell.lexa should equal(expectedLexa +- 0.00000001)
  }

  test("LexA is updated correctly after a break is repaired") {
    val l0 = LexaSteadyState / 2.0
    var cell = Cell(time = 0, lexa = l0, break = false)

    for deltaTime <- Seq(0, 1, 1.5, 3, 3.5) do
      cell = Cell.updateCell(
        reaction =
          NextReaction(reaction = EnumReactions.None, deltaTime = deltaTime),
        cell = cell
      )
      val expectedLexa =
        (l0 - LexaSteadyState) * math.exp(
          -LexaDecay * cell.time
        ) + LexaSteadyState
      cell.lexa should equal(expectedLexa +- 0.00000001)
  }
}
