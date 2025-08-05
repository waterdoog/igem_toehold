import matplotlib.pyplot as plt
import numpy as np

# Define nucleotides and their (x,y) positions manually based on the structure in the image

def plot_nucleotide(ax, x, y, base, number=None, color='black', circle_color='white'):
    """Plot a nucleotide as a circle with text inside, and optional number annotation."""
    circle = plt.Circle((x, y), 0.3, edgecolor='black', facecolor=circle_color, linewidth=1.5)
    ax.add_artist(circle)
    ax.text(x, y, base, fontsize=10, ha='center', va='center', color=color, fontweight='bold')
    if number is not None:
        ax.text(x, y-0.5, str(number), fontsize=8, ha='center', va='center', color='black')

def draw_line(ax, x1, y1, x2, y2, lw=1.5):
    """Draw a line between two points representing bonds."""
    ax.plot([x1, x2], [y1, y2], color='black', linewidth=lw)

# Nucleotide sequences from the image
# Conserved sequence loop (positions 40 to 60)
conserved_loop_bases = ['A','U','U','U','U','A','A','A','A','G','C','G','A']
conserved_loop_pos = list(range(40, 53))

# Top loop (positions 50-60 around circle at top-left)
top_loop_bases = ['A', 'G', 'A', 'G', 'G', 'A', 'G','A','C','A']
top_loop_pos = list(range(50,60))

# Stem vertical (positions 30-40 and 60-70) - All N's except annotated C and G at 40/60
stem_bases_top = ['A','G','C','G','A','U','G','G','N','N','N']  # partial example, will correct below
stem_bases = ['N']*20  # mostly 'N's from 30 to 70
stem_pos = list(range(30, 70+1))  # 30 to 70 inclusive

# Big circular loop (positions 10-30 and 70-90)
big_loop_bases_left = ['N']*21  # 10-30 (21 bases)
big_loop_bases_right = ['A','A','C','C','U','G','G','G','C','G','A','G','C','G','A','A','A','A','A']  # 70+ to 90
big_loop_pos_left = list(range(10,31))
big_loop_pos_right = list(range(80,91))

# Because sequence details are partial, we'll replicate according to the image exactly.

# Now, manually specifying positions as x,y to replicate the topology

fig, ax = plt.subplots(figsize=(10, 12))
ax.set_aspect('equal')
ax.axis('off')

# Coordinates for conserved loop: roughly circle at (0,9), radius 2
center_conserved = (0, 9)
radius_conserved = 2

# Plot conserved loop nucleotides (positions 40-60) - 13 bases on a small circle arc
theta_conserved = np.linspace(np.pi*0.25, np.pi*1.75, 13)  # arc angles from ~45deg to 315 deg

for i, (base, pos) in enumerate(zip(conserved_loop_bases, conserved_loop_pos)):
    x = center_conserved[0] + radius_conserved * np.cos(theta_conserved[i])
    y = center_conserved[1] + radius_conserved * np.sin(theta_conserved[i])
    plot_nucleotide(ax, x, y, base, number=pos, circle_color='white')

# Top small loop around positions 50-60 (outer circle at top-left)
center_toploop = (-3.5, 11)
radius_toploop = 1.3
theta_toploop = np.linspace(np.pi*0.2, np.pi*1.8, 10)

base_toploop = ['A','G','A','G','G','A','G','A','C','A']
pos_toploop = list(range(50,60))

for i, (base, pos) in enumerate(zip(base_toploop, pos_toploop)):
    x = center_toploop[0] + radius_toploop * np.cos(theta_toploop[i])
    y = center_toploop[1] + radius_toploop * np.sin(theta_toploop[i])
    plot_nucleotide(ax, x, y, base, number=pos)

# Connect top loop nucleotides inside conserved
for i in range(len(base_toploop)-1):
    x1 = center_toploop[0] + radius_toploop * np.cos(theta_toploop[i])
    y1 = center_toploop[1] + radius_toploop * np.sin(theta_toploop[i])
    x2 = center_toploop[0] + radius_toploop * np.cos(theta_toploop[i+1])
    y2 = center_toploop[1] + radius_toploop * np.sin(theta_toploop[i+1])
    draw_line(ax, x1, y1, x2, y2)

# Connect conserved loop nucleotides
for i in range(len(conserved_loop_bases)-1):
    x1 = center_conserved[0] + radius_conserved * np.cos(theta_conserved[i])
    y1 = center_conserved[1] + radius_conserved * np.sin(theta_conserved[i])
    x2 = center_conserved[0] + radius_conserved * np.cos(theta_conserved[i+1])
    y2 = center_conserved[1] + radius_conserved * np.sin(theta_conserved[i+1])
    draw_line(ax, x1, y1, x2, y2)

# Connect top loop to conserved loop (position 50)
x_top_end = center_toploop[0] + radius_toploop * np.cos(theta_toploop[-1])
y_top_end = center_toploop[1] + radius_toploop * np.sin(theta_toploop[-1])
x_cons_start = center_conserved[0] + radius_conserved * np.cos(theta_conserved[0])
y_cons_start = center_conserved[1] + radius_conserved * np.sin(theta_conserved[0])
draw_line(ax, x_top_end, y_top_end, x_cons_start, y_cons_start)

# Plot vertical stem (positions 30 to 40 and 60 to 70) with N's except some annotated bases
# From 30 to 40, mostly N's except at 40: A, C, G
# From 60 to 70, mostly N's
stem_x = 0
stem_y_positions = np.arange(2, 9.5, 0.5)  # 30 to 69 mapped linearly

# Specify bases for stem 30 to 70
stem_bases_full = ['N']*41
# positions:
# at 40: A, C, G (3 bases) mapped at indices: 10, 11, 12 (40,41,42 approx)
stem_bases_full[10] = 'A'  # position 40
stem_bases_full[11] = 'U'  #41 (not shown in image but inferred)
stem_bases_full[12] = 'A'  #42
stem_bases_full[13] = 'A'  #43
stem_bases_full[14] = 'G'  #44
stem_bases_full[15] = 'U'  #45
stem_bases_full[16] = 'A'  #46
stem_bases_full[17] = 'C'  #47
stem_bases_full[18] = 'G'  #48

# For simplicity, let's plot N from 30 to 70 with a few annotations from image

pos_stem = list(range(30, 71))
y_stem = np.linspace(2, 8, len(pos_stem))

for i, p in enumerate(pos_stem):
    base = stem_bases_full[i] if i < len(stem_bases_full) else 'N'
    plot_nucleotide(ax, stem_x, y_stem[i], base, number=p)
    if i > 0:
        draw_line(ax, stem_x, y_stem[i-1], stem_x, y_stem[i])

# Large loop at bottom (positions from 10 to 29 on the left arc, 70 to 90 on the right arc)
# We'll draw two semi-circle arcs, left arc (positions 10-29), right arc (positions 70-90)

# Left large loop
center_left_loop = (0, 1)
radius_left_loop = 7
theta_left_loop = np.linspace(np.pi/1.7, np.pi*1.3, 20)  # large arc from about 110deg to 230deg

pos_leftloop = list(range(10,30))
for i, p in enumerate(pos_leftloop):
    x = center_left_loop[0] + radius_left_loop * np.cos(theta_left_loop[i])
    y = center_left_loop[1] + radius_left_loop * np.sin(theta_left_loop[i])
    # All 'N's on left large loop except last few (24-29) marked as 'N' in the image too
    plot_nucleotide(ax, x, y, 'N', number=p)

# Right large loop
center_right_loop = (0, 1)
radius_right_loop = 7
theta_right_loop = np.linspace(-np.pi/3, np.pi/1.7, 21)  # from -60deg to 110deg

pos_rightloop = list(range(70, 91))
right_loop_bases = ['A', 'A', 'C', 'C', 'U', 'G', 'G', 'G', 'C', 'G', 'A', 'G', 'C', 'G', 'A', 'A', 'A', 'A', 'A', 'A', 'A']

for i, (p, base) in enumerate(zip(pos_rightloop, right_loop_bases)):
    x = center_right_loop[0] + radius_right_loop * np.cos(theta_right_loop[i])
    y = center_right_loop[1] + radius_right_loop * np.sin(theta_right_loop[i])
    plot_nucleotide(ax, x, y, base, number=p)

# Connect big loop nucleotides on each side
for i in range(len(pos_leftloop)-1):
    x1 = center_left_loop[0] + radius_left_loop * np.cos(theta_left_loop[i])
    y1 = center_left_loop[1] + radius_left_loop * np.sin(theta_left_loop[i])
    x2 = center_left_loop[0] + radius_left_loop * np.cos(theta_left_loop[i+1])
    y2 = center_left_loop[1] + radius_left_loop * np.sin(theta_left_loop[i+1])
    draw_line(ax, x1, y1, x2, y2)

for i in range(len(pos_rightloop)-1):
    x1 = center_right_loop[0] + radius_right_loop * np.cos(theta_right_loop[i])
    y1 = center_right_loop[1] + radius_right_loop * np.sin(theta_right_loop[i])
    x2 = center_right_loop[0] + radius_right_loop * np.cos(theta_right_loop[i+1])
    y2 = center_right_loop[1] + radius_right_loop * np.sin(theta_right_loop[i+1])
    draw_line(ax, x1, y1, x2, y2)

# Connect big loop left to stem (position 30)
x_left_end = center_left_loop[0] + radius_left_loop * np.cos(theta_left_loop[-1])
y_left_end = center_left_loop[1] + radius_left_loop * np.sin(theta_left_loop[-1])
x_stem_start = stem_x
y_stem_start = y_stem[0]
draw_line(ax, x_left_end, y_left_end, x_stem_start, y_stem_start)

# Connect big loop right to stem (position 70)
x_right_start = center_right_loop[0] + radius_right_loop * np.cos(theta_right_loop[0])
y_right_start = center_right_loop[1] + radius_right_loop * np.sin(theta_right_loop[0])
x_stem_end = stem_x
y_stem_end = y_stem[-1]
draw_line(ax, x_right_start, y_right_start, x_stem_end, y_stem_end)

# Annotations text
ax.annotate('Conserved\nsequence',
            xy=(0, 9.5), xytext=(2.5, 11.5),
            arrowprops=dict(arrowstyle='<->', lw=2),
            fontsize=12, ha='center')

ax.text(-5, 5, 'Toehold', fontsize=14, fontweight='bold', ha='center', va='center')
ax.text(5, 5, 'Universal\nlinker', fontsize=14, fontweight='bold', ha='center', va='center')
ax.text(0, -3, 'Series B', fontsize=16, fontweight='bold', ha='center', va='center')

# Set limits to nicely frame the structure
ax.set_xlim(-10, 10)
ax.set_ylim(-5, 13)

plt.show()