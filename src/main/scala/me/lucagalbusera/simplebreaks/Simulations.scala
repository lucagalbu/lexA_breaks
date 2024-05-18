package me.lucagalbusera.simplebreaks

import java.io.{File, PrintWriter}

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
