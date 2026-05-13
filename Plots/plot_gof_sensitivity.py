# Sensitivity of IC-Cox goodness-of-fit across submodel sizes (p = 8,000)
# Data source: Table 10 (d = 50); d = 15, 30, 60 from additional computation
# Methods shown use best-performing imputation strategy per baseline method
import matplotlib.pyplot as plt
import numpy as np
    
d_values = [15, 30, 50, 60]

data = {
    'Log-rank SIS': {
        'LogLik': [-357.263, -331.806, -295.865, -285.005],
        'ChiSq': [138.830, 189.743, 261.626, 283.345],
        'AIC': [744.526, 721.613, 689.730, 686.011],
        'BIC': [798.940, 826.814, 867.483, 896.413]
    },
    'ADD-SIS': {
        'LogLik': [-387.517, -369.296, -350.239, -327.011],
        'ChiSq': [78.321, 114.764, 152.878, 199.334],
        'AIC': [805.034, 796.592, 794.478, 768.022],
        'BIC': [859.449, 901.793, 964.976, 974.796]
    },
    'MV-SIS (M3)': {
        'LogLik': [-386.409, -362.292, -329.998, -322.471],
        'ChiSq': [80.538, 128.771, 193.359, 208.415],
        'AIC': [804.817, 786.585, 759.997, 762.941],
        'BIC': [862.859, 899.041, 941.378, 976.971]
    },
    'DC-SIS (M1)': {
        'LogLik': [-387.240, -360.189, -341.954, -328.795],
        'ChiSq': [78.876, 132.977, 169.447, 195.766],
        'AIC': [806.480, 782.379, 785.909, 777.590],
        'BIC': [864.522, 894.835, 970.918, 995.247]
    },
    'KS-Filter (M2)': {
        'LogLik': [-391.735, -363.075, -343.406, -329.948],
        'ChiSq': [69.885, 127.205, 166.544, 193.460],
        'AIC': [807.470, 782.150, 782.811, 773.895],
        'BIC': [851.002, 883.724, 956.937, 980.670]
    }
}

colors = {
    'Log-rank SIS': '#c44e52',
    'ADD-SIS': '#4c72b0',
    'MV-SIS (M3)': '#55a868',
    'DC-SIS (M1)': '#dd8452',
    'KS-Filter (M2)': '#8172b3'
}
markers = {
    'Log-rank SIS': 'o',
    'ADD-SIS': 'v',
    'MV-SIS (M3)': '*',
    'DC-SIS (M1)': 's',
    'KS-Filter (M2)': '^'
}

metrics = ['LogLik', 'ChiSq', 'AIC', 'BIC']
y_labels = ['Log-likelihood', 'Likelihood Ratio $\chi^2$', 'AIC', 'BIC']

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams['xtick.major.width'] = 1.8
plt.rcParams['ytick.major.width'] = 1.8

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
axes = axes.flatten()

lines = []
labels = []

for i, metric in enumerate(metrics):
    ax = axes[i]
    for method in data.keys():
        lw = 3.0 if method == 'Log-rank SIS' else 2.0
        ms = 9 if method == 'Log-rank SIS' else 8
        zorder = 10 if method == 'Log-rank SIS' else 5

        line, = ax.plot(d_values, data[method][metric],
                        marker=markers[method],
                        color=colors[method],
                        linewidth=lw,
                        markersize=ms,
                        zorder=zorder)

        if i == 0:
            lines.append(line)
            labels.append(method)

    ax.set_xticks(d_values)
    ax.set_xlabel('Submodel Size $d$', fontsize=16)
    ax.set_ylabel(y_labels[i], fontsize=16)

    ax.tick_params(axis='both', which='major', labelsize=14, direction='out', length=6)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')

fig.legend(lines, labels,
           loc='lower center',
           bbox_to_anchor=(0.5, 0.0),
           ncol=3,
           prop={'weight': 'bold', 'size': 14},
           frameon=False,
           columnspacing=2.5,
           handletextpad=0.5)

plt.tight_layout(rect=[0, 0.13, 1, 1], w_pad=0.8, h_pad=0.8)

plt.savefig('GoF_Sensitivity_Analysis_Bold_Compact.pdf', format='pdf', bbox_inches='tight')
plt.savefig('GoF_Sensitivity_Analysis_Bold_Compact.png', dpi=600, bbox_inches='tight')

plt.show()