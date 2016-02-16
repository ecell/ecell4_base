import sys
import collections

try:
    from IPython.core.display import clear_output
    have_ipython = True
except ImportError:
    have_ipython = False

class ProgressBar:

    def __init__(self, width=30, slug='#', space='-',
                 bar_template=' [{bar}]  {info}', head=None):
        self.width = width
        self.slug = slug
        self.space = space
        self.bar_template = bar_template
        self.head = head

        # self.markers = '|/-\\'
        # self.__updates = 0

        self.__last = (0.0, 0.0)

        self.update(0.0)

        if have_ipython:
            self.animate = self.animate_ipython
            self.flush = self.flush_ipython
        else:
            self.animate = self.animate_noipython
            self.flush = self.flush_noipython

    def format_eta(self, eta):
        eta = int(eta)
        seconds = eta % 60
        eta //= 60
        minutes = eta % 60
        eta //= 60
        hours = eta % 24
        eta //= 24
        if eta > 0:
            days = eta
            return '{:d} {:02d}:{:02d}:{:02d}'.format(days, hours, minutes, seconds)
        else:
            return '{:02d}:{:02d}:{:02d}'.format(hours, minutes, seconds)

    def update(self, progress, elapsed=None):
        if progress > 1.0:
            progress = 1.0
        elif progress < 0.0:
            progress = 0.0

        nmax = self.width - 2
        n = int(round(progress * nmax))

        if self.head is not None and (0 < n < nmax):
            bar = '{}{}{}'.format(self.slug * (n - 1), self.head, self.space * (nmax - n))
        else:
            bar = '{}{}'.format(self.slug * n, self.space * (nmax - n))

        # marker += '  {}'.format(self.markers[self.__updates])
        # self.__updates = (self.__updates + 1) % len(self.markers)

        info = '  {:>5.1f}%'.format(progress * 100)
        if elapsed is not None:
            info += '  Elapsed:  ' + self.format_eta(elapsed)
            # info += '  Elapsed:  ' + "{}".format(elapsed)
            if progress > self.__last[0]:
                # speed = elapsed / progress
                speed = (elapsed - self.__last[1]) / (progress - self.__last[0])
                info += ' ETA:  ' + self.format_eta(speed * (1.0 - progress))
                # info += ' ETA:  ' + "{}".format(speed * (1.0 - progress))
            self.__last = (progress, elapsed)

        items = {'bar': bar, 'info': info}
        self.progressbar = self.bar_template.format(**items)

    def flush_ipython(self):
        try:
            clear_output()
        except Exception:
            # terminal IPython has no clear_output
            pass
        print('\r', end='')
        sys.stdout.flush()

    def flush_noipython(self):
        print('\r', end='')
        sys.stdout.flush()

    def animate_ipython(self, *args, **kwargs):
        self.update(*args, **kwargs)

        try:
            clear_output()
        except Exception:
            # terminal IPython has no clear_output
            pass
        print('\r {}'.format(self.progressbar), end='')
        sys.stdout.flush()

    def animate_noipython(self, *args, **kwargs):
        self.update(*args, **kwargs)
        print('\r{}'.format(self.progressbar), end='')
        sys.stdout.flush()

class ProgressBarSimulatorWrapper:
    """A wrapper class to show a progress bar for running a simulation
    """

    def __init__(self, sim, timeout=10, flush=False, **kwargs):
        """Constructor.

        Parameters
        ----------
        sim : Simulator
            A wrapped Simulator object
        timeout : float, optional
            An interval to update the progress bar. Given as seconds.
            Default is 10.
        flush : bool, optional
            Clear the output at finishing a simulation.
            Default is False.

        See Also
        --------
        ProgressBar

        """
        self.__sim = sim
        self.__timeout = timeout
        self.__flush = flush
        self.__kwargs = kwargs

    def run(self, duration, obs):
        """Run the simulation.

        Parameters
        ----------
        duration : Real
            a duration for running a simulation.
                A simulation is expected to be stopped at t() + duration.
        observers : list of Obeservers, optional
            observers

        """
        from ecell4.core import TimeoutObserver

        timeout = TimeoutObserver(self.__timeout)
        if isinstance(obs, collections.Iterable):
            obs = tuple(obs) + (timeout, )
        else:
            obs = (obs, timeout)
        p = ProgressBar(**self.__kwargs)
        p.animate(0.0)
        tstart = self.__sim.t()
        upto = tstart + duration
        while self.__sim.t() < upto:
            self.__sim.run(upto - self.__sim.t(), obs)
            p.animate((self.__sim.t() - tstart) / duration, timeout.accumulation())
        if self.__flush:
            p.flush()
        else:
            print()

    def __getattr__(self, key):
        return getattr(self.__sim, key)

progressbar = ProgressBarSimulatorWrapper


if __name__ == "__main__":
    import time

    # print(have_ipython)

    p = ProgressBar()
    for i in range(1001):
        p.animate(i * 0.001, i * 0.03)
        time.sleep(0.03)
    print()
