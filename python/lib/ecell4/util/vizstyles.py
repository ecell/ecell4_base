import copy


def periodic_color_scale(color_list):
    class PeriodicColorScale:
        """
        Color generator
        """

        COLORS = color_list

        def __init__(self, config=None):
            """
            Initialize a color scale

            Parameters
            ----------
            config : dict, default {}
                Dict for configure default colors. Its values are colors unique
                to each key. Colors included in config will never be used.

            """
            self.__config = config or {}
            self.__buffer = copy.copy(self.COLORS)

            for color in self.__config.values():
                if color in self.__buffer:
                    self.__buffer.remove(color)
            if len(self.__buffer) == 0:
                self.__buffer = copy.copy(self.COLORS)

        def get_color(self, name):
            """
            Get color unique to the recieved name

            Parameters
            ----------
            name : string
                This method returns one color unique to this parameter.

            """
            if self.__config.get(name) is None:
                self.__config[name] = self.__buffer.pop(0)
                if len(self.__buffer) == 0:
                    self.__buffer = copy.copy(self.COLORS)

            return self.__config[name]

        def get_config(self):
            """Get an instance of dic as the config of colors."""
            return self.__config
    return PeriodicColorScale

elegans_color_scale = periodic_color_scale(
        ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#e31a1c", "#8dd3c7",
         "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69",
         "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"])
matplotlib_color_scale = periodic_color_scale(
        ["#0000ff", "#008800", "#ff0000", "#ff00ff", "#ffff00", "#00ffff",
         "#000000"])
default_color_scale = elegans_color_scale

attractive_mpl_color_scale = periodic_color_scale(
        [(0.2980392156862745, 0.4470588235294118, 0.6901960784313725), (0.3333333333333333, 0.6588235294117647, 0.40784313725490196), (0.7686274509803922, 0.3058823529411765, 0.3215686274509804), (0.5058823529411764, 0.4470588235294118, 0.6980392156862745), (0.8, 0.7254901960784313, 0.4549019607843137), (0.39215686274509803, 0.7098039215686275, 0.803921568627451)])  # originally from the seaborn color palette
