import requests
import os
from io import BytesIO
from PySide6.QtGui import QPixmap

API_KEY = os.getenv("OW_API_KEY")

def get_weather(city):
    get_weather_url = "https://api.openweathermap.org/data/2.5/weather"
    params = {
        "q": city,
        "appid": API_KEY,
    }
    response = requests.get(url = get_weather_url, params = params)
    response.raise_for_status()
    data = response.json()
    weather_id = data["weather"][0]["id"]
    weather = data["weather"][0]["description"]
    temp_in_c = round(data["main"]["temp"] - 273.15)
    return weather_id, weather, temp_in_c

icon_dict = {
    2: "11d",
    3: "09d",
    50: "10d",
    51: "13d",
    52: "09d",
    53: "09d",
    6: "13d",
    7: "50d",
    800: "01d",
    801: "02d",
    802: "03d",
    803: "04d",
    804: "04d"

}

def get_icon(id):
    weather_id_str = str(id)

    for key in icon_dict:
        if weather_id_str.startswith(str(key)):
           response = requests.get(f"https://openweathermap.org/img/wn/{icon_dict[key]}@2x.png")
           image_data = BytesIO(response.content)
           pixmap = QPixmap()
           pixmap.loadFromData(image_data.read())
           return pixmap

