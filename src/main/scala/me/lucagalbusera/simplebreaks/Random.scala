package me.lucagalbusera.simplebreaks

import scala.util.Random

trait RandomGeneratorTrait:
  def getUnifRandom(): Double
  def getExpRandom(rate: Double): Double

object RandomGenerator extends RandomGeneratorTrait:
  val randomGenerator = new Random(seed = 1234)

  def getUnifRandom() = randomGenerator.nextDouble()

  def getExpRandom(rate: Double) =
    math.log(1 - randomGenerator.nextDouble()) / (-rate)
