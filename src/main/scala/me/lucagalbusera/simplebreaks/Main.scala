package me.lucagalbusera.simplebreaks

@main def main(): Unit =
  Simulations.constitutiveExpressers(filename =
    "results/constitutiveExpressers"
  )

  Simulations.checkRepairRates(
    filename = "results/breakRepairRate"
  )

  Simulations.checkLexA(filename = "results/lexaDynamics")
