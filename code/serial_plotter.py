from serial import Serial
import numpy, time, fiddle, sys

class PowerTrace:
    def __init__(self, axes, label, t0, t, power):
        self.t0 = t0
        t -= t0
        color = next(axes._get_lines.prop_cycler)['color']
        self.line = [axes.plot([t], [power],'.', color=color)[0],
                     axes.plot([t], [power], '-',label=label, lw=3, color=color)[0]]
        axes.set_xlim(0,50)
        axes.set_ylim(-90,0)



    def add_point(self, time, power, width):
        x, y = self.line[0].get_data()
        x = numpy.append(x, time-self.t0)
        y = numpy.append(y, power)
        self.line[0].set_data(x, y)
        filter = numpy.ones(width)/width
        if width < y.size:
            self.line[1].set_data(x[width-1:], numpy.log(numpy.convolve(numpy.exp(y), filter, 'valid')))
        else:
            self.line[1].set_data([x[0],x[-1]],[numpy.average(y)]*2)


class SerialPlotter:
    def __init__(self, port):
        self.s = Serial(port=port, timeout=0.0, baudrate=9600)
        self.traces = {}
        self.t0 = time.time()

    def __call__(self, fig, width, initial=False):
        if initial:
            ax = fig.add_subplot(111)
            ax.set_xlabel('time since start (s)')
            ax.set_ylabel('power (dBm)')
        else:
            ax, = fig.get_axes()
        lines = self.s.readlines()
        for line in lines:
            line = line.decode('utf8').split()
            if len(line) != 4:
                continue
            try:
                ap = {x[0]:x[1] for x in [word.split('=') for word in line]}
            except:
                continue

            chname = 'ch ' + ap['channel'] + ": " + ap['ssid']
            power = int(ap['rssi'])
            if chname in self.traces:
                self.traces[chname].add_point(time.time(), power, width)
            else:
                self.traces[chname] = PowerTrace(ax, chname, self.t0, time.time(), power)
                ax.legend()





if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("Usage: %s ARGON_SERIAL_PORT" % sys.argv[0])
        exit(1)
    plotter = SerialPlotter(sys.argv[1])
    fiddle.fiddle(plotter, [('filter width', 1, 50, 10,1)], update_interval=0.1)
