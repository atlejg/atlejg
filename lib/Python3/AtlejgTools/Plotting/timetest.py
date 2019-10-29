import random
from datetime import datetime, timedelta

today = datetime.now()

dates = [today + timedelta(days=i) for i in range(10)]
values = [random.randint(1, 20) for i in range(10)]
plot_date(date2num(dates), values, linestyle='-')
