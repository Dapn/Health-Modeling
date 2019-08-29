
# Give the chart file a name.
png(file = "TR_Population_Distribution_no_trt.jpg")

TR <- sim_no_trt$TR # Data to plot
matplot(x=0:(length(sim_no_trt$TR[,1])-1),TR, type = c("o"),pch=20,col = 1:4,xlab = "Cycles", ylab = "Proportion Population") #plot
legend("topright", legend = c("H","S1","S2","D"), col=1:4, pch=20,) # Legend

# Save the file.
dev.off()


# Give the chart file a name.
png(file = "TR_Population_Distribution_trt.jpg")

TR <- sim_trt$TR # Data to plot
matplot(x=0:(length(sim_trt$TR[,1])-1),TR, type = c("o"),pch=20,col = 1:4,xlab = "Cycles", ylab = "Proportion Population") #plot
legend("topright", legend = c("H","S1","S2","D"), col=1:4, pch=20,) # Legend

# Save the file.
dev.off()

