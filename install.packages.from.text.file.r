########## INSTALL PACKAGES FROM TEXT FILE ##############

installPackagesTxt = function(file, nomultiarch = TRUE, ask = F, repos = 1:10){
  
  pkgs = readLines(file)
  pkgs = pkgs[-c(grep("^#",pkgs))]
  pkgs = pkgs[pkgs!=""]
  
  print(pkgs)
  new.pkgs <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  
  print(new.pkgs)
  setRepositories(ind = repos)

  
  if(length(new.pkgs)){
	for(i in 1:length(new.pkgs)){
   
		print(paste("package ",i,": ",new.pkgs[i],sep=""))
			if(nomultiarch){
				install.packages(new.pkgs[i], INSTALL_opts = c("--no-multiarch"))
			} else {
				install.packages(new.pkgs[i])
			}
     
		}
	}
}