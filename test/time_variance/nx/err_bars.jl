using Plots

# Sample data
x = 1:5
y = [2.0, 3.0, 4.0, 5.0, 6.0]

# Asymmetric errors: lower and upper errors
yerr_lower = [0.5, 0.4, 0.3, 0.2, 0.1]
yerr_upper = [0.2, 0.3, 0.4, 0.5, 0.6]

# Create the plot with asymmetric error bars
plot(x, y, yerror = (yerr_lower, yerr_upper), label = "Data with asymmetric error bars")

# Add titles and labels
title!("Asymmetric Error Bars Example")
xlabel!("X-axis")
ylabel!("Y-axis")

# Display the plot
show()
