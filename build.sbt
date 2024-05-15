val scala3Version = "3.4.1"

lazy val root = project
  .in(file("."))
  .settings(
    name := "basic breaks simulations",
    version := "0.1.0",
    scalaVersion := scala3Version,
    libraryDependencies += "org.scalamock" %% "scalamock" % "6.0.0" % Test,
    libraryDependencies += "org.scalatest" %% "scalatest" % "3.2.18" % "test",
    scalacOptions += "-Yrangepos"
  )
