package me.lucagalbusera.simplebreaks

import java.io.{File, PrintWriter}
import scala.collection.mutable.ListBuffer

object Simulations:
  def constitutiveExpressers(filename: String): Unit =
    val file = new File(filename)
    val bw = new PrintWriter(file)

    bw.write(headerString() + "\n")

    val rates = Rates(
      mrnaCreation = 10,
      mrnaDecay = 1,
      breakRepair = 0,
      proteinFolding = 1,
      lexaUnbinding = LexaUnbinding,
      proteinCreation = ProteinCreation
    )
    var cell = Cell(repressed = false, lexa = 0)
    var reaction = NextReaction(reaction = EnumReactions.None, deltaTime = 0)

    while cell.time < 50 do
      bw.write(
        infoString(cell = cell, rates = rates, reaction = reaction) + "\n"
      )

      reaction = NextReaction.computeNext(
        rates = rates,
        cell = cell,
        randomGenerator = RandomGenerator
      )
      cell = Cell.updateCell(reaction = reaction, cell = cell)
    bw.close()

  def checkRepairRates(filename: String): Unit =
    val file = new File(filename)
    val bw = new PrintWriter(file)

    bw.write("iter" + headerString() + "\n")

    for repairRate <- Seq(0.1, 0.5, 1) do
      println(s"Repair rate: $repairRate")

      val rates = Rates.createFromSlope(
        mrnaDecay = 0.25,
        breakRepair = repairRate,
        proteinFolding = 20
      )

      for iter <- 0 to 500 do
        var cell = Cell()
        var reaction =
          NextReaction(reaction = EnumReactions.None, deltaTime = 0)
        var repairedBreaks = 0

        while reaction.reaction != EnumReactions.breakRepaired do
          reaction = NextReaction.computeNext(
            rates = rates,
            cell = cell,
            randomGenerator = RandomGenerator
          )
          cell = Cell.updateCell(reaction = reaction, cell = cell)
          bw.write(
            s"$iter" +
              infoString(cell = cell, rates = rates, reaction = reaction) +
              "\n"
          )
          if reaction.reaction == EnumReactions.breakRepaired then
            repairedBreaks += 1

    bw.close()

  /** Simulates a cell without breaks and reports only the binding/unbinding of
    * LexA. In this way, the binding/unbinding rates can be checked.
    *
    * @param filename
    *   Name of a file where to write the results.
    */
  def checkLexaRates(filename: String): Unit =
    val file = new File(filename)
    val bw = new PrintWriter(file)

    bw.write(headerString() + "\n")

    val rates = Rates(
      mrnaCreation = 2,
      mrnaDecay = 0.25,
      breakRepair = 0,
      proteinFolding = 1.0 / 7.0,
      lexaUnbinding = LexaUnbinding,
      proteinCreation = 20
    )

    var cell = Cell(break = false)
    var reaction =
      NextReaction(reaction = EnumReactions.None, deltaTime = 0)
    var numLexaReactions = 0

    while numLexaReactions < 50000 do
      reaction = NextReaction.computeNext(
        rates = rates,
        cell = cell,
        randomGenerator = RandomGenerator
      )
      cell = Cell.updateCell(reaction = reaction, cell = cell)
      if reaction.reaction == EnumReactions.lexaBinding || reaction.reaction == EnumReactions.lexaUnbinding
      then
        numLexaReactions += 1
        bw.write(
          infoString(cell = cell, rates = rates, reaction = reaction) +
            "\n"
        )

    bw.close()

  def checkLexA(filename: String): Unit =
    val file = new File(filename)
    val bw = new PrintWriter(file)

    bw.write(headerString() + "\n")

    val rates = Rates.createFromSlope(
      mrnaDecay = 0.25,
      breakRepair = 0.5,
      proteinFolding = 20
    )

    var cell = Cell()
    var reaction =
      NextReaction(reaction = EnumReactions.None, deltaTime = 0)

    while cell.time < 30 do
      bw.write(
        infoString(cell = cell, rates = rates, reaction = reaction) + "\n"
      )
      reaction = NextReaction.computeNext(
        rates = rates,
        cell = cell,
        randomGenerator = RandomGenerator
      )
      cell = Cell.updateCell(reaction = reaction, cell = cell)

    bw.close()

  /** Simulate various breaks and computes their heights and widths. The height
    * is defined as the highest point of the gfp production (protein *
    * protein_folding), while the width is the full width at half maximum.
    *
    * @param filename
    *   Name of a file where to write the results. The full simulation is saved
    *   to filename_full, while the peaks stats are saved to filename_stats.
    * @param rates
    *   Constat rates to use for the simulations.
    * @param numBreaks
    *   Number of breaks to simulate.
    * @param threshold
    *   Height at which the signal protein*rf is considered to be in a peak.
    * @param printFullResults
    *   Should all the simulation be printed (filename_full), or only the peaks
    *   stats (filename_stats)? Useful when we want to simulate a large number
    *   of breaks.
    */
  def breakDurationsFull(
      filename: String,
      rates: Rates,
      numBreaks: Int,
      threshold: Double,
      printFullResults: Boolean
  ): Unit =
    val fileFull = new File(filename + "_full")
    val bwFull = new PrintWriter(fileFull)
    if printFullResults then bwFull.write("peakNum;" + headerString() + "\n")

    val fileStats = new File(filename + "_stats")
    val bwStats = new PrintWriter(fileStats)
    bwStats.write("peakNum;rr;dm;rf;height;duration\n")

    var foundBreaks = 0

    while foundBreaks < numBreaks do
      var cell = Cell()
      var reaction =
        NextReaction(reaction = EnumReactions.None, deltaTime = 0)

      var measures = ListBuffer[Cell](cell)
      while cell.break || cell.protein > 0 do
        reaction = NextReaction.computeNext(
          rates = rates,
          cell = cell,
          randomGenerator = RandomGenerator
        )
        cell = Cell.updateCell(reaction = reaction, cell = cell)
        measures += cell

      val height = measures.map(cell => cell.protein).max
      if height * rates.proteinFolding / CellVolume > threshold then
        println(s"Found peaks: $foundBreaks")
        if printFullResults then
          measures.foreach(cell =>
            bwFull.write(
              s"$foundBreaks;" +
                infoString(cell = cell, rates = rates, reaction = reaction) +
                "\n"
            )
          )
        val height_index = measures
          .map(cell => cell.protein)
          .indexOf(height)
        val start = measures
          .take(height_index + 1)
          .findLast(cell => cell.protein <= height / 2.0)
        val end = measures
          .drop(height_index)
          .findLast(cell => cell.protein >= height / 2.0)
        val duration =
          (end.get.time - start.get.time) / (2.0 * math.sqrt(
            2.0 * math.log(2.0)
          ))

        bwStats.write(
          s"$foundBreaks;${rates.breakRepair};${rates.mrnaDecay};${rates.proteinFolding};" +
            s"${height * rates.proteinFolding / CellVolume};$duration\n"
        )
        foundBreaks += 1

    if printFullResults then bwFull.close()
    bwStats.close()

  private def headerString(): String =
    "rm;dm;rp;rf;rr;ld;lss;loff;time;reaction;lexa;break;repressed;mrna;protein"

  private def infoString(
      cell: Cell,
      rates: Rates,
      reaction: NextReaction
  ): String =
    s"${rates.mrnaCreation};${rates.mrnaDecay};${rates.proteinCreation};" +
      s"${rates.proteinFolding};${rates.breakRepair};" +
      s"${LexaDecay};${LexaSteadyState};${rates.lexaUnbinding};" +
      s"${cell.time};${reaction.reaction};${cell.lexa};${cell.break};${cell.repressed};" +
      s"${cell.mrna};${cell.protein}"
