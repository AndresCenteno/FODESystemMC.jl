using DataFrames
using CSV
using Plots

df = CSV.read("test/time_variance/nx/nx_alpha_T_1D_heat_results.csv")

# four plots, plot xlabel nx, for many alpha
# Your code to fill the DataFrame

# Assuming you've filled the DataFrame `df` with appropriate values

# Group the DataFrame by alpha
grouped_df = groupby(df, :alpha)

# Plot for each alpha
plot()
for g in grouped_df
    alpha_val = first(g[:alpha])
    nx_vals = g[:nx]
    mean_t_vals = g[:mean_t]
    plot!(nx_vals, mean_t_vals, label = "Alpha = $alpha_val", marker = :circle)
end

# Customize plot
xlabel!("nx")
ylabel!("Mean Time")
title!("Mean Time vs nx for different alpha values")
display(plot!)