#Add warningsFileName as a global variable so it doesn't have to be passed around as much
utils::globalVariables("warningsFileName")
#Environment to add .RPPA.fit.model to so we don't have to make a call to assign to the GlobalEnvironment
RPPASPACE_Temp_Env <- new.env();