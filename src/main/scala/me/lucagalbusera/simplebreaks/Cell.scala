package me.lucagalbusera.simplebreaks

case class Cell(
    time: Double = 0,
    lexa: Double = 200,
    break: Boolean = true,
    repressed: Boolean = true,
    mrna: Int = 0,
    protein: Int = 0
)

object Cell:
  def updateCell(reaction: NextReaction, cell: Cell): Cell =
    val newTime = cell.time + reaction.deltaTime
    val newLexa = updateLexa(
      nextTime = reaction.deltaTime,
      break = cell.break,
      currentLexa = cell.lexa,
      lexaDecay = LexaDecay,
      lexaSteadyState = LexaSteadyState
    )

    reaction.reaction match
      case EnumReactions.None =>
        Cell(
          mrna = cell.mrna,
          protein = cell.protein,
          time = newTime,
          break = cell.break,
          repressed = cell.repressed,
          lexa = newLexa
        )
      case EnumReactions.mrnaCreation =>
        Cell(
          mrna = cell.mrna + 1,
          protein = cell.protein,
          time = newTime,
          break = cell.break,
          repressed = cell.repressed,
          lexa = newLexa
        )
      case EnumReactions.mrnaDecay =>
        Cell(
          mrna = cell.mrna - 1,
          protein = cell.protein,
          time = newTime,
          break = cell.break,
          repressed = cell.repressed,
          lexa = newLexa
        )
      case EnumReactions.proteinCreation =>
        Cell(
          mrna = cell.mrna,
          protein = cell.protein + 1,
          time = newTime,
          break = cell.break,
          repressed = cell.repressed,
          lexa = newLexa
        )
      case EnumReactions.proteinFolded =>
        Cell(
          mrna = cell.mrna,
          protein = cell.protein - 1,
          time = newTime,
          break = cell.break,
          repressed = cell.repressed,
          lexa = newLexa
        )
      case EnumReactions.breakRepaired =>
        Cell(
          mrna = cell.mrna,
          protein = cell.protein,
          time = newTime,
          break = false,
          repressed = cell.repressed,
          lexa = newLexa
        )
      case EnumReactions.lexaBinding =>
        Cell(
          mrna = cell.mrna,
          protein = cell.protein,
          time = newTime,
          break = cell.break,
          repressed = true,
          lexa = newLexa
        )
      case EnumReactions.lexaUnbinding =>
        Cell(
          mrna = cell.mrna,
          protein = cell.protein,
          time = newTime,
          break = cell.break,
          repressed = false,
          lexa = newLexa
        )

  private def updateLexa(
      nextTime: Double,
      break: Boolean,
      currentLexa: Double,
      lexaDecay: Double,
      lexaSteadyState: Double
  ) =
    if (break) then currentLexa * math.exp(-nextTime * lexaDecay)
    else
      (currentLexa - lexaSteadyState) * math.exp(
        -lexaDecay * nextTime
      ) + lexaSteadyState
