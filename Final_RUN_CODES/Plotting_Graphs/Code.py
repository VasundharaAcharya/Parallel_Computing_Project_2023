# importing pandas library
import pandas as pd
# import matplotlib library
import matplotlib.pyplot as plt
import textwrap
  
# creating dataframe
df = pd.read_csv("MPI_CODE_VERSUS_SERIAL_CODE.csv")

#print(df.head())
# plotting graph
ax=df.plot(x="Code", y=["ClockCycles"], kind="bar")
wrapped_labels = [textwrap.fill(label, 8) for label in df["Code"]]
# set font size of x-axis and y-axis labels
ax.set_xlabel("Code", fontsize=20,fontweight='bold')
ax.set_ylabel("Clockcycles", fontsize=20,fontweight='bold')
ax.set_xlim([-0.5, 3.5])
# set font size of x-tick and y-tick labels
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_xticklabels(wrapped_labels, rotation=0)
# set font size of legend
ax.legend(fontsize=12, title_fontsize=12, bbox_to_anchor=(1.05, 1))

# show the plot
plt.show()
