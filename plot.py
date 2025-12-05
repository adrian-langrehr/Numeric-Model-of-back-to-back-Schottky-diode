# Load a Result object and make a simple plot. For more sophisticated plots result.plot() is helpful.
import dill
import matplotlib.pyplot as plt

with open('example.dill', 'rb') as f:
    result = dill.load(f)

print(result.fit_report())

V = result.userkws['V']

plt.figure(figsize=(5,4))
plt.plot(result.userkws['V'], result.data*1e9, 'o', label='Data')
plt.plot(result.userkws['V'], result.best_fit*1e9, '-', label='Fit')

dely = result.eval_uncertainty(x=V, sigma=3)
plt.fill_between(V, (result.best_fit-dely)*1e9, (result.best_fit+dely)*1e9, color='#888888', alpha=0.5, label='Confidence Interval')
plt.legend()
plt.xlabel('Voltage (V)')
plt.ylabel(r'Current (nA)')
plt.grid()
plt.tight_layout()
plt.show()