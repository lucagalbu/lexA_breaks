val scala3Version = "3.4.1"

lazy val root = project
  .in(file("."))
  .settings(
    name := "basic breaks simulations",
    version := "0.1.0",
    scalaVersion := scala3Version,
    libraryDependencies += "org.scalameta" %% "munit" % "1.0.0-M12" % Test
  )
