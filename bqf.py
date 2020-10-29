from sage.all import *
from matplotlib import pyplot as plt
from matplotlib import animation
from math import sqrt, sin, cos, acos

STEPS = 50
PAUSE = 10

def get_lattice(cur_bqf):
    coeffs = cur_bqf.coefficients()

    base = sqrt(float(coeffs[0]))
    side = sqrt(float(coeffs[2]))

    phi = acos(float(coeffs[1])/sqrt(float(coeffs[0])*float(coeffs[2])))
    
    x1 = [0.,side*cos(phi)]
    y1 = [0.,side*sin(phi)]
    
    x2 = [0., base]
    y2 = [0.,0.]

    return x1, y1, x2, y2, phi

#deforms the identity matrix smoothly to the desired matrix
def slow_transform(R, cur_step, steps=STEPS):
    s = float(cur_step)
    steps = float(steps)
    
    b = 0.*(1.-(s/steps)) + float(R[0,1])*(s/steps)
    c = 0.*(1.-(s/steps)) + float(R[1,0])*(s/steps)
    d = 1.*(1.-(s/steps)) + float(R[1,1])*(s/steps)
    a = (b*c + 1)/d

    #print(R)
    #print(R[0][1])
    #print(s/steps)
    #print(Matrix([[a,b],[c,d]]))
    
    return Matrix([[a,b],[c,d]])

def update(num, bqf, R, line1, line2, cur_phi, steps=STEPS):
    if num <= steps:
        #deforming the bqf slightly
        cur_bqf = bqf(slow_transform(R, num, steps=STEPS))
        #print(cur_bqf)
        x1, y1, x2, y2, phi = get_lattice(cur_bqf)

        #print("Step {} out of {}".format(num+1, steps))

        line1.set_data(x1,y1)
        line2.set_data(x2,y2)
        cur_phi.set_text('Phi: {:.2f}'.format(phi))

    else:
        pass

bqf = QuadraticForm(RR, 2, [1,2,5])

bqf.compute_definiteness()
if not bqf.is_positive_definite():
    raise Exception("Gotta be pos. def. yo")

disc = bqf.disc()

#R is the matrix which reduces the bqf
G, R = bqf.reduced_binary_form1()

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title("Reduction of the BQF {:.2f}x^2 + {:.2f}xy + {:.2f}y^2".format(float(bqf.coefficients()[0]),float(bqf.coefficients()[1]),float(bqf.coefficients()[2])))
plt.xlim([-.1,max(bqf.coefficients())])
plt.ylim([-.1,max(bqf.coefficients())])

x1, y1, x2, y2, phi = get_lattice(bqf)

line1, = ax.plot(x1,y1, color='b')
line2, = ax.plot(x2,y2, color='b')
cur_phi = ax.text(max(bqf.coefficients())-1., max(bqf.coefficients())-1., 'Phi: {:.2f}'.format(phi), fontsize=15)

ani = animation.FuncAnimation(fig, update, STEPS+PAUSE, fargs=(bqf, R, line1, line2, cur_phi,), blit=False)
ani.save('reducing_bqf.gif', writer='imagemagick')

print("Original BQF:")
print(bqf.coefficients())
print("Reduced BQF:")
print(G)
