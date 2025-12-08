from genice3.genice import GenIce3
from logging import getLogger, basicConfig, INFO

logger = getLogger("test_genice3api1")
basicConfig(level=INFO)
genice = GenIce3()
logger.info("Reactive properties:")
logger.info(f"     All: {genice.list_all_reactive_properties().keys()}")
logger.info(f"  Public: {genice.list_public_reactive_properties().keys()}")
logger.info("Settabe reactive properties:")
logger.info(f"     All: {genice.list_settable_reactive_properties().keys()}")
logger.info(f"  Public: {genice.list_public_settable_reactive_properties().keys()}")
