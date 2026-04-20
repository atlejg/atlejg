#################### EXAMPLE DEMO <ex_demo.py> ###############################
'''A simple interactive demo to illustrate the use of IPython's Demo class.'''

print('Hello, welcome to an interactive IPython demo.')

# The mark below defines a block boundary, which is a point where IPython will
# stop execution and return to the interactive prompt. The dashes are actually
# optional and used only as a visual aid to clearly separate blocks while
# editing the demo code.
# <demo> stop

x = 1
y = 2

# <demo> stop

# the mark below makes this block as silent
# <demo> silent

print('This is a silent block, which gets executed but not printed.')

# <demo> stop
# <demo> auto
print('This is an automatic block.')
print('It is executed without asking for confirmation, but printed.')
z = x + y

print('z =', x)

# <demo> stop
# This is just another normal block.
print('z is now:', z)

print('bye!')
################### END EXAMPLE DEMO <ex_demo.py> ############################
