import org.scalatest.funspec._
import org.scalatest.matchers.should._
import org.scalatest.funsuite.AnyFunSuite
import me.lucagalbusera.simplebreaks.RandomGenerator

class ExponentialDistributionTest extends AnyFunSuite with Matchers {
  test("summary statistics of generated numbers are correct within 5%") {
    val numSamples = 10000

    for rate <- Seq(0.5, 3) do {
      val generatedValues =
        (1 to numSamples).map(_ => RandomGenerator.getExpRandom(rate = rate))

      // Test mean
      val mean = generatedValues.sum.toDouble / numSamples
      val mean_expected = 1.0 / rate
      mean should equal(mean_expected +- mean_expected * 0.05)

      // Test variance
      val variance =
        generatedValues.map(x => Math.pow(x - mean, 2)).sum / numSamples
      val variance_expected = 1.0 / (rate * rate)
      variance should equal(variance_expected +- variance_expected * 0.05)
    }
  }
}
