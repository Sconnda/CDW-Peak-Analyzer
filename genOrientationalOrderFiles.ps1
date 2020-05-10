$alpha = 0.06
$size = 768,1024,1536,2048

For ($i = 0; $i -lt $size.Length; $i++) {
	For ($j = 0; $j -lt $alpha.Length; $j++) {
		python3 LatticeRelaxation.py $size[$i] $alpha[$j]
		python3 OrientationalOrder.py SimulatedData
		cd SimulatedData

		$name = -join ($alpha[$j],"_",$size[$i],"_RelaxationAnimation.avi")
		Rename-Item -Path "RelaxationAnimation.avi" -NewName $name
		$name = -join ($alpha[$j],"_",$size[$i],"_SimulatedData_G6.csv")
		Rename-Item -Path "SimulatedData_G6.csv" -NewName $name

		$path = -join ($alpha[$j],"_",$size[$i],"_RelaxationAnimation.avi")
		$destination = -join ("G6 Files/",$alpha[$j],"_",$size[$i],"_RelaxationAnimation.avi")
		Move-Item -Path $path -Destination $destination
		$path = -join ($alpha[$j],"_",$size[$i],"_SimulatedData_G6.csv")
		$destination = -join ("G6 Files/",$alpha[$j],"_",$size[$i],"_SimulatedData_G6.csv")
		Move-Item -Path $path -Destination $destination

		Remove-Item SimulatedData*

		cd ..
	}
}