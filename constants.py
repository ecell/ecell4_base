class EventType(object):
    def __init__(self, value, symbol):
        self.value = value
        self.symbol = symbol

    def __int__(self):
        return self.value

    def __repr__(self):
        return "EventType(%s, %s)" % (repr(self.value), repr(self.symbol))

EventType.SINGLE_REACTION = EventType(10, "SINGLE_REACTION")
EventType.SINGLE_ESCAPE   = EventType(11, "SINGLE_ESCAPE")
EventType.COM_ESCAPE      = EventType(11, "COM_ESCAPE")
EventType.IV_EVENT        = EventType(12, "IV_EVENT")
EventType.IV_ESCAPE       = EventType(13, "IV_ESCAPE")
EventType.IV_REACTION     = EventType(14, "IV_REACTION")
EventType.IV_INTERACTION  = EventType(15, "IV_INTERACTION")
EventType.BURST           = EventType(16, "BURST")
EventType.MULTI_ESCAPE    = EventType(17, "MULTI_ESCAPE")
EventType.MULTI_REACTION  = EventType(18, "MULTI_REACTION")
